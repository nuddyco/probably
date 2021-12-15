(defpackage :probably
  (:use :common-lisp :alexandria :serapeum :let-plus)
  (:shadow #:trace))
(in-package :probably)

(defmacro λ (args &body body) `(lambda ,args ,@body))

;;; TODO
;;; - Separate get-random-number and logprob
;;; - Metropolis-Hastings
;;; - Gamma distribution

;;;; Basic types

(deftype model ()
  "Probabilistic model.
   A Lisp function that samples random variables."
  'function)

(deftype choicemap ()
  "Chosen values of named random variables."
  'hash-table)

(deftype logprob ()
  "Logarithm of a probability value.
   Used for the sake of floating point accuracy."
  '(real * #.(log 1)))

(deftype probability ()
  "Probability value in [0,1]."
  '(real 0 1))

(deftype addr ()
  "Address of a random value in a trace.
   Can have any concrete type."
  t)

(defunit ◯
  "Special value for 'missing' values of unbound random variables.")

(defun present? (value)
  "True if VALUE is not the missing value ◯."
  (not (eq value ◯)))

;;;; Traces

(defstruct trace
  "Recording of one execution of a model.
   Records the random choices made and their joint probability."
  (model       (required-argument) :type model)
  (args        nil                 :type list)
  (choices     (dict)              :type choicemap)
  (constraints (dict)              :type choicemap)
  (score       (log 1.0d0)         :type logprob))

(defvar-unbound *trace*
  "The trace that is currently being recorded.")

(defvar *observations* (dict)
  "Fixed random choices provided from outside the model.")

(defun bound? (addr &optional (trace *trace*))
  "Is the random variable with ADDR already bound to a value?"
  (hash-table-includes? (trace-choices trace) addr))

(defun binding (addr &optional (trace *trace*))
  (@ (trace-choices trace) addr))

(defun constrained? (addr &optional (trace *trace*))
  "Is the random variable with ADDR constrained to a specific value?"
  (hash-table-includes? (trace-constraints trace) addr))

(defun constraint (addr &optional (trace *trace*))
  (@ (trace-constraints trace) addr))

(-> hash-table-includes (hash-table t) (member t nil))
(defun hash-table-includes? (hash-table key)
  "Does HASH-TABLE include KEY?"
  (nth-value 1 (@ hash-table key)))

;;;; ADTs

;;;;; Model

(-> generate (model &key (:args t) (:constraints choicemap)) (values t probability trace))
(defun generate (model &key args (constraints (dict)))
  (let ((*trace* (make-trace :model model :args args :constraints constraints)))
    (values (apply model args) (exp (trace-score *trace*)) *trace*)))

;;;;; Trace

(-> logpdf (trace) number)
(defun logpdf (trace)
  (trace-score trace))

(-> choices (trace) choicemap)
(-> update (trace choicemap &optional t) trace)

;;;; Probability distributions

(-> bernoulli (number addr &optional t) t)
(defun bernoulli (p addr &optional (outcome ◯))
  (probably addr outcome
            (op (bernoulli-sample p)) (op (bernoulli-log-likelihood p _))))

(defun bernoulli-sample (p)
  "Return one sample from a Bernoulli(P) distribution."
  (< (random 1.0) p))

(defun bernoulli-log-likelihood (p outcome)
  "Return the Bernoulli(P) log-likelihood of OUTCOME."
  (log (if outcome p (- 1 p))))

(-> probably (addr t function function) t)
(defun probably (addr observation sample-fn log-likelihood-fn)
  "Bind a random variable to its value and record probability score."
  (let ((value (cond ((present? observation) observation)
                     ((constrained? addr)    (constraint addr))
                     (t                      (funcall sample-fn)))))
    (when (bound? addr)
      ;; Maybe we should allow compatible values? Not sure yet. Play it safe.
      (with-simple-restart (continue "Continue with value ~s" value)
        (error "Address ~s already bound to ~s when binding to ~s")))
    (bind! addr value (funcall log-likelihood-fn value))
    value))

(-> bind! (addr t logprob) t)
(defun bind! (addr value log-likelihood)
  "Bind ADDRESS to VALUE accounting for LOG-LIKELIHOOD."
  ;; XXX Gen registers prefixes of addresses. Why?
  (incf (trace-score *trace*) log-likelihood)
  (setf (@ (trace-choices *trace*) addr) value))

;;;; Importance sampling
;;;
;;; Bayesian Data Analysis (BDA3), section 10.4.

#|
(-> importance-sampling (model choicemap model (integer 0 *)))
(defun importance-sampling (model observations proposal-model n)
  (loop for i from 0 below n 
        for proposed-trace = (simulate proposal-model)
        for all-choices = (merge observations (trace-choices proposed-trace))
        for trace = (generate model all-choices)
        for weight = (- (logpdf trace) (logpdf proposed-trace))
        collecting trace into traces
        collecting (log weight) into log-weights
        finally (let ((ℓ (logsumexp log-weights)))
                  (return (mapcar (λ (trace log-weight)
                                    (list trace (- log-weight ℓ)))
                                  traces log-weights)))))

 (return (mapcar (λ (trace weight)
                                  (list trace (/ weight (sum weights)))))

        collecting trace into traces
        collecting (- (logpdf trace) (logpdf proposed-trace)) into weights
        finally (return (mapcar (λ (trace weight)
                                  (list trace (- weight (logsumexp weights))))))
        summing (logpdf proposed-trace) into sum
        collect (list trace weight)))

(defun logsumexp (seq)
  ;; Equations 1.6 and 1.7 of Gen thesis
  (let ((max (extremum seq #'>)))
    (+ max (log (loop for x in seq summing (exp (- x max)))))))
|#

;;;; Metropolis-Hastings

(defun metropolis-hastings (model n &key args (warmup 100) (constraints (dict)))
  "Return N representative samples from MODEL."
  (subseq (loop with sample = (lambda ()
                                (generate model :args args :constraints constraints))
                with remaining = (+ n warmup)
                with p₀ = (nth-value 1 (funcall sample))
                for (value p trace) = (multiple-value-list (funcall sample))
                for accept-threshold = (/ p p₀)
                when (> accept-threshold (random 1.0))
                  collect (prog1 trace
                            (setf p₀ p)
                            (decf remaining))
                while (plusp remaining))
          warmup))

;;;; Testing

(defun burglary ()
  "Burglary model as a probabilistic program.
   :BURGLARY - did a burglary occur?
   :DISABLED - was the alarm disabled?
   :ALARM    - did the alarm sound?
   :CALLS    - did the phone ring?"
  (let* ((burglary (bernoulli 0.01 :burglary))
         (disabled (and burglary
                        (bernoulli 0.1 :disabled)))
         (alarm (and (not disabled)
                     (bernoulli (if burglary 0.94 0.01) :alarm)))
         (calls (bernoulli (if alarm 0.70 0.05) :calls)))
    calls))

(defun investigate-calls (&optional (samples 1000000))
  (let* ((traces (metropolis-hastings #'burglary samples
                                      :constraints (dict :calls t)))
         (burglaries (count-if (op (binding :burglary _)) traces)))
    (format t "~&Out of ~A calls there were ~A burglaries (~,2f%)~%"
            samples burglaries (* 100.0 (/ burglaries samples)))))

(defun investigate-burglaries (&optional (samples 1000000))
  (let* ((traces (metropolis-hastings #'burglary samples
                                      :constraints (dict :burglary t)))
         (calls (count-if (op (binding :calls _)) traces)))
    (format t "~&Out of ~A burglaries there were ~A calls (~,2f%)~%"
            samples calls (* 100.0 (/ calls samples)))))

