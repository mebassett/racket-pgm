#lang typed/racket
(require "model.rkt" "random-var.rkt" "factors.rkt")


; step 1 - let's define some discrete random variables.
; the built-in constructer isn't great, use make-DiscreteRandomVar as it
; does some housekeeping.

(define difficulty (make-DiscreteRandomVar 'difficulty '("hard" "easy")))
(define grade (make-DiscreteRandomVar 'grade '("a" "b" "f")))
(define intelligence (make-DiscreteRandomVar 'iq '("high" "low")))
(define SAT (make-DiscreteRandomVar 'SAT '("high" "low")))
(define letter (make-DiscreteRandomVar 'letter '("good" "poor")))


; okay, now let's make it a directed graphical model by adding some dependencies!

(add-dependency! grade difficulty)
(add-dependency! grade intelligence)
(add-dependency! SAT intelligence)
(add-dependency! letter grade)

; now we can call this a model.
; again, don't use the built-in constructor for DiscreteModel.
;
; make-naive-model initializes everything to use a uniform distribution across
; all labels.
;
; this whole model thing is really just to keep track of which random vars
; belong to which CPDs and such.  It could probably be moved into a custodian.
; (which smells a lot like a probalistic programming language.)
;

(define model (make-naive-model (list difficulty grade SAT letter intelligence)))

; now, we don't really want uniform distributions.  So we will updated each
; cpd with values from the book.

(update-cpd! model difficulty (make-factor-from-table* difficulty #hash((((difficulty . "hard")) . 0.4) (((difficulty . "easy")) . 0.6))))
(update-cpd! model intelligence (make-factor-from-table* intelligence #hash((((iq . "high")) . 0.3) (((iq . "low")) . 0.7))))
(update-cpd! model grade (make-factor-from-table* grade
#hash(
       (((difficulty . "easy") (iq . "low") (grade . "a")) . 0.3)
       (((grade . "b") (difficulty . "easy") (iq . "low")) . 0.4)
       (((grade . "f") (difficulty . "easy") (iq . "low")) . 0.3)
       (((difficulty . "hard") (iq . "low") (grade . "a")) . 0.05)
       (((difficulty . "hard") (grade . "b") (iq . "low")) . 0.25)
       (((difficulty . "hard") (grade . "f") (iq . "low")) . 0.7)
       (((iq . "high") (difficulty . "easy") (grade . "a")) . 0.9)
       (((iq . "high") (grade . "b") (difficulty . "easy")) . 0.08)
       (((iq . "high") (grade . "f") (difficulty . "easy")) . 0.02)
       (((iq . "high") (difficulty . "hard") (grade . "a")) . 0.5)
       (((iq . "high") (difficulty . "hard") (grade . "b")) . 0.3)
       (((iq . "high") (difficulty . "hard") (grade . "f")) . 0.2)
       )))
(update-cpd! model SAT (make-factor-from-table* SAT
#hash(
       (((SAT . "low") (iq . "low")) . 0.95)
       (((SAT . "high") (iq . "low")) . 0.05)
       (((iq . "high") (SAT . "low")) . 0.2)
       (((iq . "high") (SAT . "high")) . 0.8)
       )))
(update-cpd! model letter (make-factor-from-table* letter
    #hash(
       (((letter . "poor") (grade . "a")) . 0.1)
       (((letter . "good") (grade . "a")) . 0.9)
       (((grade . "b") (letter . "poor")) . 0.4)
       (((grade . "b") (letter . "good")) . 0.6)
       (((grade . "f") (letter . "poor")) . 0.99)
       (((grade . "f") (letter . "good")) . 0.01)
       )))

; so the data structures I'm using here are pretty wonky.
; it would probably be best if I had some notion of an instance of
; a set of variables at particular values.
; instead I am keeping the entire set of all possible values around all the time.
; this gets really weird in the source for Factors.

; Anyhow.  Onwards and upwards.

;(is-valid-model? model)

; let's hope that returns true!  That just checks to make sure all CPDs exists
; and are sane.  Let's try some inferrence now:

; what if we want to know the probability that a student will fail given that
; he did well on his SATs and it's a really hard course.
((probability-conditioned-on-evidence model (set grade) (set '(difficulty . "hard")
                                                             '(SAT . "high")))
 (set '(grade . "f")))

; should be 0.26363636363636367 or thereabouts

; you can see there that I am using a (Pairof RandomVar Symbol) to keep track of
; which variable is initiated to what.
; you can also see that the student does not have too high of a probablity of failing.

; let's try another one.
; what about the probability the course is hard
; given that he got a good letter and has good SATs?

((probability-conditioned-on-evidence model (set difficulty)
                                      (set '(SAT . "high")
                                           '(letter . "good")))
 (set '(difficulty . "hard")))

; should be 0.3209440506070644 or thereabouts.
; that's all for now!

