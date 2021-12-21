;;; probably.lisp – a toy probabilistic programming framework based on Gen.jl
;;
;; This code was written as an exercise to understand what
;; probabilistic programming is all about, particularly while looking
;; at Gen.jl (http://gen.dev/) and reading part of the associated thesis.
;;
;; Written by Luke Gorrie <luke@nuddy.co> in 2021. MIT license.

(defpackage :probably
  (:use :common-lisp :alexandria :serapeum)
  (:shadow #:trace))
(in-package :probably)

(eval-always (setf *read-default-float-format* 'double-float))

(defmacro λ (args &body body) `(lambda ,args ,@body))

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
  '(double-float * #.(log 1d0)))

(deftype probability ()
  "Probability value in [0,1]."
  '(double-float 0d0 1d0))

(deftype addr ()
  "Address of a random value in a trace.
   Can have any concrete type."
  t)

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

(def □ '□
  "Special value for 'missing' values of unbound random variables.")

(deftype □ () '(eq □))

(defun present? (value)
  "True if VALUE is not the missing value □."
  (not (eq value □)))

;;;; Core engine

;;;;; Making random choices

(-> probably (addr t function function) t)
(defun probably (addr observation sample-fn log-likelihood-fn)
  "Bind a random variable to its value and record probability score."
  (let ((value (cond ((present? observation) observation)
                     ((constrained? addr)    (constraint addr))
                     (t                      (funcall sample-fn)))))
    (assert (present? value))
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

;;;;; Running models

(-> generate (model &key (:args t) (:constraints choicemap)) (values t logprob trace))
(defun generate (model &key args (constraints (dict)))
  "Generate a trace by simulating MODEL.
   CONSTRAINTS optionally maps random variable addresses to fixed values."
  (let ((*trace* (make-trace :model model :args args :constraints constraints)))
    (values (apply model args) (trace-score *trace*) *trace*)))


(defun metropolis-hastings (model n &key args (warmup 100) (constraints (dict)))
  "Return N representative samples from MODEL."
  (sort (subseq (loop for seq from 0
                      with sample = (lambda ()
                                      (generate model :args args :constraints constraints))
                      with collected = (- warmup)
                      with ref = (nth-value 1 (funcall sample))
                      for (value logprob trace) = (multiple-value-list (funcall sample))
                      for α = (if (> logprob ref)
                                  1d0
                                  (exp (max (log double-float-epsilon) (- logprob ref))))
                      for r = (random 1.0d0)
                      when (<= r α)
                        collect (prog1 trace
                                  (setf ref logprob)
                                  (incf collected))
                while (< collected n))
                warmup)
        #'> :key #'trace-score))

;;;; Querying and updating traces

(defun bound? (addr &optional (trace *trace*))
  "Is the random variable with ADDR already bound to a value?"
  (hash-table-includes? (trace-choices trace) addr))

(defun bindings (&optional (trace *trace*))
  (hash-table-plist (trace-choices trace)))

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

;;;; Probability distributions

;;;;; Bernoulli

(-> bernoulli (number addr &optional t) t)
(defun bernoulli (p addr &optional (outcome □))
  (probably addr outcome
            (op (bernoulli-sample p))
            (op (bernoulli-log-likelihood p _))))

(defun bernoulli-sample (p)
  "Return one sample from a Bernoulli(P) distribution."
  (< (random 1.0d0) p))

(defun bernoulli-log-likelihood (p outcome)
  "Return the Bernoulli(P) likelihood of OUTCOME."
  (log (if outcome p (- 1 p))))

;;;;; Normal

;;(-> normal (addr real real &optional (or real □)) t)
(defun normal (addr μ σ &optional (outcome □))
  (probably addr outcome
            (op (normal-sample μ σ))
            (op (normal-log-likelihood μ σ _))))

(defun normal-sample (μ σ)
  "Sample one value from Normal(μ,σ)."
  (+ μ (* σ (unit-normal-sample))))

(defun unit-normal-sample ()
  "Sample one value from Normal(0,1)."
  (loop for u1 = (random 1.0d0)
        for u2 = (random 2.0d0)
        when (> u1 single-float-epsilon)
          do (return (* (sqrt (* -2 (log u1)))
                        (cos (* pi 2 u2))))))

(defun normal-log-likelihood (μ σ x)
  "Return the likelihood of Normal(μ,σ) at X."
  (let ((z (/ (- x μ) σ)))
    (- 0
       (log σ)
       (/ (+ (* z z) (log (* pi 2))) 2))))

;;;;; Known

(defun known (addr value)
  "Known VALUE at ADDRESS with probability 1."
  (probably addr value
            (op value)
            (op (if (eq _ value) (log 1d0) 0d0))))

;;;; Examples

;;;;; Burglary model

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

;;;;; Linear regression

(defun center (xs)
  (loop with μ = (normal :μ (mean xs) 1)
        for i from 0
        for x in xs
        do (normal `(:x ,i) μ 1 x)))

(defun linear-model (xs)
  ;; Model parameters with broad priors
  (let ((slope     (normal :slope     0 10))
        (intercept (normal :intercept 0 10)))
    ;; Parameter (:X <index>) for each value in XS.
    (loop for i from 0
          for x in xs
          for prediction = (+ intercept (* slope i))
          ;; Predict points to be near the line
          ;; described by the model parameters.
          do (normal `(:x ,i) prediction 0.1 x))))

(defun report (traces)
  (format t "~&intercept   slope   log-probability")
  (loop for trace in (sort (copy-list traces) #'< :key #'trace-score)
        do (format t "~& ~8,3f  ~6,3f  ~10,5f"
                   (binding :intercept trace)
                   (binding :slope trace)
                   (trace-score trace))))

;;;;; Polynomial regression

(defun polynomial-model (xs)
  (loop ;; Keep adding powers while a coin-toss comes up heads.
        with order = (known :order (loop while (bernoulli 0.5 (gensym)) count t))
        ;; Distribute coefficients normally around zero.
        with coeff = (loop for i upto order collect (normal i 0 0.5))
        ;; The polynomial prediction function
        with polyf  = (λ (n) (loop for k in coeff
                                   for power from 0
                                   summing (* k (expt n power))))
        ;; Make prediction
        for n from 0
        for prediction = (funcall polyf n)
        ;; Match prediction to observation
        for x in xs
        do (normal `(:x ,n) prediction 1 x)))

(defun polyreport (traces)
  (let ((freq (frequencies traces :key (op (binding :order _)))))
    (format t "~&How many terms (n) in k₀ + k₁*x + k₂*x² + ... + kₙ*xⁿ?~%~
                 Terms: ~{~4d ~}~%~
                 Count: ~{~4d ~}~%"
            (iota 10)
            (mapcar (op (or (@ freq _) 0)) (iota 10)))))

(defun polydemo ()
  ;; might take a minute or so...
  (let* ((data (loop for i from 1 to 5
                     collect (+ 1.0d0 (* i -1) (* 0.5 (expt i 3)))))
         (traces (metropolis-hastings (op (polynomial-model data)) 100)))
    (polyreport traces)))

    
