#lang typed/racket
(require math/flonum)

(struct RandomVar ([name : Symbol]
                   [depends-on : (Listof RandomVar)]
                   [influences : (Listof RandomVar)])
  #:mutable)
(struct DiscreteRandomVar RandomVar ([labels : (Listof String)]))

(define (add-dependency! [var : RandomVar]
                         [dependency : RandomVar]) : Void
  (set-RandomVar-influences! dependency (cons var (RandomVar-influences dependency)))
  (set-RandomVar-depends-on! var (cons dependency (RandomVar-depends-on var))))

(define (make-DiscreteRandomVar [name : Symbol]
                                [labels : (Listof String)]) : DiscreteRandomVar
  (DiscreteRandomVar name '() '() labels))


(define-type TableCPDIndex (Pairof Symbol String))
(define-type TableCPD (HashTable (Setof TableCPDIndex)
                                 Float))

(define (make-TableCPD-labels [vars : (Listof DiscreteRandomVar)]) : (Listof (Listof TableCPDIndex))
  (apply cartesian-product (map (λ ([var : DiscreteRandomVar])
                                  (map (λ ([label : String])
                                         (cons (RandomVar-name var) label))
                                       (DiscreteRandomVar-labels var)))
                                vars)))
(define (get-variable-labels [var : DiscreteRandomVar]) : (Setof (Setof TableCPDIndex))
  (list->set
   (for/list : (Listof (Setof TableCPDIndex))
     ([i (make-TableCPD-labels (filter DiscreteRandomVar?
                                       (cons var (RandomVar-depends-on var))))])
     (list->set i))))


(define (generate-uniform-table-cpd [var : DiscreteRandomVar]) : TableCPD
  (for/hash : TableCPD ([k (make-TableCPD-labels (filter DiscreteRandomVar?
                                                         (cons var (RandomVar-depends-on var))))])
    (values (list->set k)
            (exact->inexact (/ 1 (length (DiscreteRandomVar-labels var)))))))

(struct DiscreteModel ([vars : (Listof DiscreteRandomVar)]
                       [cpds : (HashTable DiscreteRandomVar TableCPD)]) #:mutable)

(define (is-valid-cpd? [var : DiscreteRandomVar]
                       [cpd : TableCPD]) : Boolean
  (define (all-cells-exists?) : Boolean
    (equal? (list->set (hash-keys cpd))
            (list->set (for/list : (Listof (Setof (Pair Symbol String)))
                         ([l (make-TableCPD-labels (filter DiscreteRandomVar?
                                                           (cons var (RandomVar-depends-on var))))])
                         (list->set l)))))
  (define (rows-sum-to-one?) : Boolean
    (define rows (make-TableCPD-labels (filter DiscreteRandomVar? (RandomVar-depends-on var))))
    (define (sum-cell-row [pre-cellkeys : (Listof TableCPDIndex)]) : Float
      (flsum (map (λ ([cellkey : (Listof TableCPDIndex)]) (hash-ref cpd (list->set cellkey)))
                  (map (λ ([col : String])
                         (cons (cons (RandomVar-name var) col) pre-cellkeys))
                       (DiscreteRandomVar-labels var)))))
    (andmap (λ ([row : (Listof TableCPDIndex)]) (= (sum-cell-row row) 1)) rows))
  (and (all-cells-exists?)
       (rows-sum-to-one?)))
;
(define (is-valid-model? [model : DiscreteModel]) : Boolean
  (define (all-cpds-exist?)
    (for/and : Boolean ([v (DiscreteModel-vars model)]) (not (false? (hash-ref (DiscreteModel-cpds model) v #f)))))
  (define (all-cpds-valid?)
    (for/and : Boolean ([(var cpd) (DiscreteModel-cpds model)]) (is-valid-cpd? var cpd)))
  (and (all-cpds-exist?)
       (all-cpds-valid?)))

(define (make-naive-model [vars : (Listof DiscreteRandomVar)]) : DiscreteModel
  (define cpds : (HashTable DiscreteRandomVar TableCPD) (make-hash))
  (for ([var vars])
    (hash-set! cpds var (generate-uniform-table-cpd var)))
  (DiscreteModel vars cpds))

(define (update-cpd! [model : DiscreteModel]
                     [var : DiscreteRandomVar]
                     [pre-cpd : (HashTable (Listof TableCPDIndex) Float)]) : Void
  (define cpd
    (for/hash : TableCPD ([([k : (Listof TableCPDIndex)]
                            [v : Float])          pre-cpd])
      (values (list->set k) v)))
  (unless (member var (DiscreteModel-vars model))
    (error "Random Variable ~a not in Model ~a" var model))
  (unless (is-valid-cpd? var cpd)
    (error "TableCPD ~a is not valid for random variable ~a.  Please check its dependencies and that it sums to 1." cpd var))
  (hash-set! (DiscreteModel-cpds model) var cpd))

(define-struct/exec DiscreteFactor
  ([scope : (Setof (Setof TableCPDIndex))]
   [func : (-> (Setof TableCPDIndex) Float)])
  [(λ (self labels)
     (unless (set-member? (DiscreteFactor-scope self) labels)
        (error "The variables ~a are outside of the scope of Factor ~a" labels self (DiscreteFactor-scope self)))
     ((DiscreteFactor-func self) labels)) : (-> DiscreteFactor (Setof TableCPDIndex) Float)])

(define (associated-factor [var : DiscreteRandomVar]
                             [model : DiscreteModel]) : DiscreteFactor
    (define cpd (hash-ref (DiscreteModel-cpds model) var))
    (DiscreteFactor (list->set (hash-keys cpd)) (λ (labels) (hash-ref cpd labels))))

(define (var-in-scope? [var : DiscreteRandomVar] [factor : DiscreteFactor]) : Boolean
  (for/and : Boolean
    ([i (DiscreteFactor-scope (associated-factor grade model))]
     [j (get-variable-labels difficulty)])
    (subset? j i)))

(define (factor-marginalization [factor : DiscreteFactor]
                                [var : DiscreteRandomVar]
                                [model : DiscreteModel]) : DiscreteFactor
  (define var-scope (list->set (hash-keys (hash-ref (DiscreteModel-cpds model) var)))) 
  (unless (var-in-scope? var factor)
    (error "Variable ~a is not in the scope of Factor ~a" var factor))
  (DiscreteFactor (for/set : (Setof (Setof TableCPDIndex))
                    ([i (DiscreteFactor-scope (associated-factor grade model))]
                     [j (get-variable-labels difficulty)])
                    (set-subtract i j))
                  (λ (factor-labels)
                    (flsum (map (λ ([var-labels : (Setof TableCPDIndex)]) : Float
                                  (factor (set-union factor-labels var-labels)))
                                (set->list (get-variable-labels var)))))))

(define (get-scope-symbols [scope : (Setof (Setof TableCPDIndex))]) : (Setof Symbol)
    (list->set (set-map (first (set->list scope)) (λ ([x : TableCPDIndex]) (car x)))))

(define (labels-within-scope [scope-symbols : (Setof Symbol)]
                             [labels : (Setof TableCPDIndex)]) : (Setof TableCPDIndex)

      (list->set
       (filter (λ ([x : TableCPDIndex]) (set-member? scope-symbols (car x)))
               (set->list labels))))

(define (product-factor [f1 : DiscreteFactor]
                        [f2 : DiscreteFactor]) : DiscreteFactor
  (define f1-symbols (get-scope-symbols (DiscreteFactor-scope f1)))
  (define f2-symbols (get-scope-symbols (DiscreteFactor-scope f2)))
  (define joint-scope 
    (for/set : (Setof (Setof TableCPDIndex))
      ([i (DiscreteFactor-scope f1)]
       [j (DiscreteFactor-scope f2)])
      (set-union i j)))
  (define (factor-func [labels : (Setof TableCPDIndex)]) : Float
    (* (f1 (labels-within-scope f1-symbols labels))
       (f2 (labels-within-scope f2-symbols labels))))
  (DiscreteFactor joint-scope factor-func))
  


(define difficulty (make-DiscreteRandomVar 'difficulty '("hard" "easy")))
(define grade (make-DiscreteRandomVar 'grade '("a" "b" "f")))
(define SAT (make-DiscreteRandomVar 'SAT '("high" "low")))
(define letter (make-DiscreteRandomVar 'letter '("good" "poor")))
(define intelligence (make-DiscreteRandomVar 'iq '("high" "low")))

(add-dependency! grade difficulty)
(add-dependency! grade intelligence)
(add-dependency! SAT intelligence)
(add-dependency! letter grade)
(define model (make-naive-model (list difficulty grade SAT letter intelligence)))

(update-cpd! model difficulty #hash((((difficulty . "hard")) . 0.4) (((difficulty . "easy")) . 0.6)))
(update-cpd! model intelligence #hash((((iq . "high")) . 0.3) (((iq . "low")) . 0.7)))
(update-cpd! model grade 
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
       ))
(update-cpd! model SAT 
#hash(
       (((SAT . "low") (iq . "low")) . 0.95)
       (((SAT . "high") (iq . "low")) . 0.05)
       (((iq . "high") (SAT . "low")) . 0.2)
       (((iq . "high") (SAT . "high")) . 0.8)
       ))
(update-cpd! model letter
    #hash(
       (((letter . "poor") (grade . "a")) . 0.1)
       (((letter . "good") (grade . "a")) . 0.9)
       (((grade . "b") (letter . "poor")) . 0.4)
       (((grade . "b") (letter . "good")) . 0.6)
       (((grade . "f") (letter . "poor")) . 0.99)
       (((grade . "f") (letter . "good")) . 0.01)
       ))




                                   