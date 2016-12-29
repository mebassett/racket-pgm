#lang typed/racket
(provide make-graph
         make-empty-graph
         Graph
         has-vertex?
         has-edge?
         vertex=? 
         add-edge!
         remove-edge!
         add-vertex!
         remove-vertex!
         rename-vertex!
         get-vertices
         get-leaves
         in-vertices
         get-neighbors
         in-neighbors
         get-edges
         in-edges
         transpose
         graph-copy)

(define-type (AdjacencyList t) 
 (HashTable t (Setof t)))

(: add-edge@ (∀ (t) (-> (AdjacencyList t) t t Void)))
(define (add-edge@ adj u v)
 (hash-update! adj
               u
               (λ ([vs : (Setof t)]) (set-add vs v))
               (λ () (ann (set) (Setof t)))))

(: remove-edge@ (∀ (t) (-> (AdjacencyList t) t t Void)))
(define (remove-edge@ adj u v)
 (hash-update! adj u (λ ([vs : (Setof t)]) (set-remove vs v))))

(: add-vertex@ (∀ (t) (-> (AdjacencyList t) t Void)))
(define (add-vertex@ adj v)
 (hash-update! adj
               v
               (λ ([vs : (Setof t)]) vs)
               (λ () (ann (set) (Setof t)))))

(: remove-vertex@ (∀ (t) (-> (AdjacencyList t) t Void)))
(define (remove-vertex@ adj v)
 (hash-remove! adj v)
 (for ([(vertex edges) (in-hash adj)])
  (hash-set! adj vertex (set-remove edges v))))

(: rename-vertex@ (∀ (t) (-> (AdjacencyList t) t t Void)))
(define (rename-vertex@ adj old new)
 (hash-set! adj new (hash-ref adj old))
 (hash-remove! adj old)
 (for ([(u vs) (in-hash adj)] #:when (set-member? vs old))
  (hash-set! adj u (set-add (set-remove vs old) new))))


(struct (t) Graph ([adj : (AdjacencyList t)]))

(: get-vertices (∀ (t) (-> (Graph t) (Listof t))))
(define (get-vertices g)
 (hash-keys (Graph-adj g)))

(: get-leaves (∀ (t) (-> (Graph t) (Listof t))))
(define (get-leaves g)
 (filter (λ ([v : t]) (null? (get-neighbors g v)))
         (get-vertices g))) 

(: in-vertices (∀ (t) (-> (Graph t) (Sequenceof t))))
(define (in-vertices g)
 (in-hash-keys (Graph-adj g)))

(: get-neighbors (∀ (t) (-> (Graph t) t (Listof t))))
(define (get-neighbors g v)
 (sequence->list (in-graph-neighbors g v)))

(: in-neighbors (∀ (t) (-> (Graph t) t (Sequenceof t))))
(define (in-neighbors g v)
 (in-graph-neighbors g v))

(: vertex=? (∀ (t) (-> (Graph t) t t Boolean)))
(define (vertex=? g u v)
 (equal? u v))

(: add-edge! (∀ (t) (-> (Graph t) t t Void)))
(define (add-edge! g u v)
 (define adj (Graph-adj g))
 (add-edge@ adj u v)
 (add-vertex@ adj v))

(: remove-edge! (∀ (t) (-> (Graph t) t t Void)))
(define (remove-edge! g u v)
 (define adj (Graph-adj g))
 (remove-edge@ adj u v))

(: add-vertex! (∀ (t) (-> (Graph t) t Void)))
(define (add-vertex! g v)
 (add-vertex@ (Graph-adj g) v))

(: remove-vertex! (∀ (t) (-> (Graph t) t Void)))
(define (remove-vertex! g v)
 (remove-vertex@ (Graph-adj g) v))

(: rename-vertex! (∀ (t) (-> (Graph t) t t Void)))
(define (rename-vertex! g old new)
 (when (member new (get-vertices g))
  (error 'rename-vertex! 
         "new vertex ~a already exists in graph ~g"
         new g))
 (rename-vertex@ (Graph-adj g) old new))

(: has-vertex? (∀ (t) (-> (Graph t) t Boolean)))
(define (has-vertex? g v)
 (member v (get-vertices g)) #t)

(: has-edge? (∀ (t) (-> (Graph t) t t Boolean)))
(define (has-edge? g u v)
 (and (has-vertex? g u) 
      (has-vertex? g v)
      ; typed/racket funkiness.  member does not return Boolean
      (not (false? (member v (get-neighbors g u))))))

(: in-edges (∀ (t) (-> (Graph t) (Sequenceof (Pairof t t)))))
(define (in-edges g)
 (in-list (get-edges g)))

(: get-edges (∀ (t) (-> (Graph t) (Listof (Pairof t t)))))
(define (get-edges g)
 (for*/list : (Listof (Pairof t t)) 
            ([u (in-vertices g)] 
            [v (in-neighbors g u)])
  (cons u v)))

(: graph-copy (∀ (t) (-> (Graph t) (Graph t))))
(define (graph-copy g)
 (struct-copy Graph g [adj (hash-copy (Graph-adj g))]))

(: transpose (∀ (t) (-> (Graph t) (Graph t))))
(define (transpose g)
 (define adj^T : (AdjacencyList t) (make-hash))
 (for ([u (in-vertices g)])
  (add-vertex@ adj^T u)
  (for ([v (in-neighbors g u)])
   (add-edge@ adj^T v u)))
 (Graph adj^T))





(: in-graph-neighbors (∀ (t) (-> (Graph t) t (Sequenceof t))))

(define (in-graph-neighbors g v)
 (in-set (hash-ref (Graph-adj g) v
                   (λ () (error 'in-vertices 
                                "vertex ~a not in graph ~a"
                                v g)))))

; types are not first-class citizens (l'o'l)
; so I couldn't think of a good way to do this as a function.
(define-syntax-rule (make-empty-graph Type)
    (Graph (ann (make-hash) (AdjacencyList Type))))

(: make-graph (∀ (t) (-> (Listof (Pairof t (Option t))) (Graph t))))
(define (make-graph es)
 (define adj : (AdjacencyList t) (make-hash))
 (for ([e : (Pairof t (Option t)) es])
  (let ([u (car e)]
        [v (cdr e)])
   (cond [v (add-edge@ adj u v)
            (add-vertex@ adj v)]
         [else (add-vertex@ adj u)])))
 (Graph adj))
