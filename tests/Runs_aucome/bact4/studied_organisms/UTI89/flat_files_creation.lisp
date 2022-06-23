;; Most files should begin by specifying that they will be interpreted
;; within the Lisp package called ecocyc.
(in-package :ecocyc)
;; Select the organism PGDB as the current PGDB
(select-organism :org-id 'UTI89)
;; Create attribute-values files without the progression pop-up
(let ((*progress-noter-enabled?* NIL))
        (create-flat-files-for-current-kb))