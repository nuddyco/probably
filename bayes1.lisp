;;; bayes1.lisp – toy probabilistic programming for Bayesian inference

(defpackage :bayes1
  (:use :common-lisp :alexandria :serapeum)
  (:nicknames :bayes))
(in-package :bayes1)

;;;; Hidden state tracking
;;;
;;; Every time a random variable is assigned a value, either by
;;; sampling at random or by observing an actual value from outside
;;; the model, the overall joint probability must be updated.

(defvar *joint-probability* nil
  "Joint probability of all the random choices that have been made.")

(defun update-probability (p)
  "Update joint probability after a random choice with likelihood P."
  (setf *joint-probability* (* p *joint-probability*)))

;;;; Bernoulli distribution
;;;
;;; Support drawing random boolean values by "tossing a biased coin."
;;; The coin has probability P of coming up heads (T) and probability
;;; 1-P of coming up tails (NIL). The parameter P is the bias of the
;;; coin, which would be 0.5 for a "fair" coin.
;;;
;;; https://en.wikipedia.org/wiki/Bernoulli_distribution

(defun bernoulli (p &optional (outcome (bernoulli-sample p)))
  "Choose a random boolean value and update joint probability.
   The value is either OUTCOME or a Bernoulli random sample.

   The joint probability is always updated to reflect the probability
   of the assigned value, regardless of how it was chosen."
  (update-probability (bernoulli-likelihood p outcome))
  outcome)

(defun bernoulli-sample (p)
  "Return T with probability P and NIL with probability (- P 1)."
  (< (random 1.0d0) p))

(defun bernoulli-likelihood (p outcome)
  "Return the likelihood of OUTCOME from Bernoulli(P)."
  (if outcome p (- 1 p)))

;;;; Probabilistic program: Coin biased estimator.
;;;
;;; OBSERVE-FLIPS is a probabilistic program – a model.

(defun observe-flips (outcomes)
  "Pick a uniform random bias and observe OUTCOMES."
  (loop ;; Choose a random bias.
        ;; FIXME: Cheating to call random directly but ~okay for Uniform?
        with bias = (random 1.0)
        ;; Register each observation from outside the model as the
        ;; outcome of a Bernoulli draw with the chosen bias.
        ;;
        ;; (If this outcome is unlikely then that will be reflected in
        ;; the updated value of *joint-probability*.)
        for outcome in outcomes
        do (bernoulli bias outcome)
        finally (return bias)))

;;;; Sampling from models
;;;
;;; Sampling from the model is simple: we run it a whole bunch of
;;; times with a mixture of observed and randomized values, and each
;;; time – for the OBSERVE-FLIPS model specifically – it provides us
;;; with a possible bias of the coin and the relative likelihood of
;;; that value with respect to the observed data.
;;;
;;; Then we can make this representative of probabilities by drawing a
;;; large number of such samples and discarding most of the less
;;; probably ones in a systematic way.
;;;
;;; This is the Metropolis-Hastings MCMC algorithm:
;;;   https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

(defun simulate (model)
  "Simulate one run of MODEL and return the joint probability."
  (let ((*joint-probability* 1.0d0))
    (values (funcall model) *joint-probability*)))

(defun sample (model n &optional (warmup 100))
  "Return N representative samples from MODEL."
  (subseq (loop with remaining = (+ n warmup)
                with p₀ = (nth-value 1 (simulate model))
                for (value p) = (multiple-value-list (simulate model))
                for accept-probability = (/ p p₀)
                when (> accept-probability (random 1.0))
                  collect (prog1 value
                            (setf p₀ p)
                            (decf remaining))
                while (plusp remaining))
          warmup))

;;;; Reporting

(defun summarize (samples &key (step 0.05) (stream t))
  "Print a histogram summarizing probable coin biases based on SAMPLES."
  ;; Super quick hack to summarize the posterior distribution.
  (loop for bucket from 0 below 1 by step
        for count = (count-if (lambda (x)
                                (and (<= bucket x)
                                     (< x (+ bucket step))))
                              samples)
        do (format stream "~&~4,2f-~4,2f: ~a~%"
                   bucket (+ bucket step)
                   (make-string (floor (* 100.0 (/ count (length samples))))
                                :initial-element #\*))))

