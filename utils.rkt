#lang typed/racket

(provide set-filter)

(: set-filter (All (T) (-> (-> T Boolean) (Setof T) (Setof T))))
(define (set-filter valid? set)
  (list->set (filter valid? (set->list set))))