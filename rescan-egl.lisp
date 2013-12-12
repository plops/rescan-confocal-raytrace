;; simulator for stephans rescanning confocal microscope use maxima to
;; solve constraints for mirror placement eventually, it should be
;; possible to investigate how image distortions are related to the
;; focal lengths and geometry of the reimaging path

;; Copyright 2013 Martin Kielhorn 
    
;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

#.(load "~/quicklisp/setup")
(eval-when (:execute :compile-toplevel :load-toplevel) 
  (ql:quickload "external-program")
  (ql:quickload :bordeaux-threads))


(eval-when (:execute :compile-toplevel :load-toplevel) 
  (setf asdf:*central-registry*
	(union 
	 (list *default-pathname-defaults* "/home/pi/raspberry-pi-egl-test/")
	 asdf:*central-registry*))
  (asdf:load-system :egl)
  (asdf:load-system :cl-opengl))


(defpackage :g (:use :cl :gl :ccl))
(in-package :g)

(ccl:use-interface-dir "raspberry-pi-vc")

(deftype vec ()
  `(simple-array double-float (3)))

(defun make-vec (&optional (x 0d0) (y 0d0) (z 0d0))
  (declare ;(type double-float x y z)
	   #+sbcl (values vec &optional))
  (make-array 3 
	      :element-type 'double-float
	      :initial-contents (mapcar #'(lambda (x) (* 1d0 x)) (list x y z))))

(defmacro v (&rest args)
  `(make-vec ,@(mapcar #'(lambda (x)
			   (* 1d0 x))
		       args)))

(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
                        `(,array (&rest indices) `(aref ,',array ,@indices)))
                      arrays)
     ,@body))

(defun v. (a b)
  (declare (type vec a b)
	   #+sbcl (values double-float &optional))
  (let ((sum 0d0))
    (declare (type double-float sum))
    (with-arrays (a b)
      (dotimes (i 3)
	(incf sum (* (aref a i) (aref b i)))))
    sum))

#+nil
(v. (v 1 2 3) (v 2 3 4))


(defmacro def-v-op (op)
  `(defun ,(alexandria:format-symbol :g "V~a" op) (a b)
     (declare (type vec a b)
	      #+sbcl (values vec &optional))
     (let ((result (v)))
       (with-arrays (result a b)
	(dotimes (i 3)
	  (setf (result i) (,op (a i) (b i)))))
       result)))

(def-v-op +)
(def-v-op -)

#+nil
(v- (v) (v 1 2 3))
#+nil
(v+ (v 2 3 4) (v 1 2 3))

(defun v* (scalar v)
  (declare (type double-float scalar)
	   (type vec v)
	   #+sbcl (values vec &optional))
  (let ((result (v)))
    (with-arrays (result v)
      (dotimes (i 3)
	(setf (result i) (* scalar (v i)))))
    result))

#+nil
(v* 2d0 (v 1 2 3))

(defun vx (a b)
  "Cross product"
  (declare (type vec a b)
	   #+sbcl (values vec &optional))
  (with-arrays (a b)
    (make-vec (- (* (a 1) (b 2))
		 (* (a 2) (b 1)))
	      (- (* (a 2) (b 0))
		 (* (a 0) (b 2)))
	      (- (* (a 0) (b 1))
		 (* (a 1) (b 0))))))
#+nil
(vx (v 1)
    (v 0 1))

(defun norm (v)
  (declare (type vec v)
	   #+sbcl (values double-float &optional))
  (let ((l2 (v. v v)))
    (declare (type double-float l2))
    (sqrt l2)))

(defun normalize (v)
  (declare (type vec v)
	   #+sbcl (values vec &optional))
  (let ((l (norm v)))
    (if (zerop l)
	(v 0 0 1)
	(v* (/ l) v))))

(deftype mat ()
  `(simple-array double-float (3 3)))

(defun m (a b c d e f g h i)
  (make-array '(3 3)
              :element-type 'double-float
              :initial-contents (list (list a b c)
				      (list d e f)
				      (list g h i))))

(defun rotation-matrix (angle vec)
  "Create matrix that rotates by ANGLE radians around the direction
 VECT. VECT must be normalized."
  (declare ;(type double-float angle)
	   (type vec vec)
	   #+sbcl (values mat &optional))
  (with-arrays (vec)
    (let* ((u (vec 0)) (v (vec 1)) (w (vec 2))
	   (c (cos (* pi (/ 180) angle)))
	   (s (sin (* pi (/ 180) angle)))
	   (1-c (- 1 c))
	   (su (* s u)) (sv (* s v)) (sw (* s w)))
      (m (+ c (* 1-c u u))   (+ (* 1-c u v) sw)   (- (* 1-c u w) sv)
	 (- (* 1-c u v) sw)  (+ c (* 1-c v v))    (+ (* 1-c v w) su)
	 (+ (* 1-c u w) sv)  (- (* 1-c v w) su)   (+ c (* 1-c w w))))))
#+nil
(rotation-matrix 90 (v 0 0 1))

(defun m* (m v)
  "Multiply matrix M with vector V."
  (declare (type mat m)
	   (type vec v)
	   #+sbcl (values vec &optional))
  (let ((res (v)))
    (with-arrays (res v m)
      (dotimes (i 3)
	(dotimes (j 3)
	  (incf (res i) (* (m i j) (v j))))))
    res))
#+nil
(m* (rotation-matrix 90 (v 0d0 0d0 1d0)) (v 1d0))


;; f1 f2 A B C D
(defparameter *maxima-batch*
  "--batch-string=load(minpack)$
f1:~A;
f2:~A;
ll:2*(f1+f2);
aa:~A;
bb:~A;
cc:~A;
dd:~A;
gg:aa*ll;
d:bb*gg;
a:cc*gg;
eq:[
	a+b+c+d+sqrt(d^2+h1^2)+f+g+sqrt(a^2+h2^2)-ll,
	a+b+c+d-gg,
	h1/c-h2/b,
	sqrt(h1^2+c^2)-f,
	sqrt(h2^2+b^2)-g,
	h1/f-h2/g,
	dd-h1/h2]; 
sol:minpack_lsquares(eq,[b,c,h1,h2,f,g],[1,1,1,1,1,1]);
/* output ll,gg,d,h1,a,h2 */
solution:[ll,gg,d,sol[1][3],a,sol[1][4]];
")

(defun get-parameters-from-maxima (&key (f1 50) (f2 100) (a .2) (b .2) (c .2) (d 1))
 (let ((s (make-string-output-stream)))
  #+sbcl (sb-ext:run-program "/usr/local/bin/maxima" 
		       (list (format nil *maxima-batch* f1 f2 a b c d)) :output s)
  #-sbcl (external-program:run "/usr/bin/maxima" 
			       (list (format nil *maxima-batch* f1 f2 a b c d)) :output s)
  (let* ((maxima-output (get-output-stream-string s))
	  (line-start (position #\Newline maxima-output 
				:from-end t :end (1- (length maxima-output))))
	  (line-end (position #\Newline maxima-output :from-end t))
	  (line (subseq maxima-output (1+ line-start) line-end))
	  (data (subseq line (+ 1 (position #\[ line)) (position #\] line)))
	  (n (1+ (count-if #'(lambda (x) (char= x #\,)) data)))
	  (num-start 0))
     (loop for i below n collect
	  (multiple-value-bind (val new-num-start) 
	      (read-from-string data nil nil :start num-start)
	    (setf num-start (1+ new-num-start))
	    val)))))
#+nil
(defparameter *geometry-parameters* (get-parameters-from-maxima :a .25 :b .1 :c .2 :d 1))

(defvar *geometry-parameters* (list 300
				    75.0
				    7.5
				    52.492676
				    15.0
				    52.492676))

(defvar *angle1* 0)
(defvar *angle2* 0)
(defvar *vertex-buf* (make-heap-ivector (* 3 2) 'single-float))
(defun store-vertex-list (vs)
  (when (< (length *vertex-buf*) (* 3 (length vs)))
    (dispose-heap-ivector *vertex-buf*)
    (make-heap-ivector (* 3 (length vs)) 'single-float))
  (let ((i 0)) 
    (loop for p in vs do 
	 (loop for c across p do
	      (setf (aref *vertex-buf* i) (coerce c 'single-float))
	      (incf i)))))

(defvar *initialized* nil)
(let ((v 0))
  (defun bla ()
    (incf v 9)
    (when (< 360 v)
      (setf v 0))
    (clear-color (* .5 (+ 1 (sin (* pi (/ v 360f0))))) 0 0 .4)
    (clear :color-buffer-bit)
    (load-identity)
    (translate 0 0 -1)
    (scale .1 .1 .1)
    (rotate v 0 0 1)
    (unless (and *geometry-parameters* 
		 (find-if-not #'numberp *geometry-parameters*))
      (destructuring-bind (ll gg d h1 a h2) *geometry-parameters*
	(declare (ignore ll))
	(let* ((alpha 10)
	       (m-sys (rotation-matrix alpha (v  0 0 1)))
	       (d1 (m* m-sys (make-vec gg)))
	       (d2 (m* m-sys (make-vec (- d) (- h1))))
	       (d3 (m* m-sys (make-vec (- (- gg (+ d a))) (+ h1 h2))))
	       (d4 (m* m-sys (make-vec (- a) (- h2))))
	       (v (v 0 0))
	       (vs (loop for e in (list (v 0 0) (v 0 10)
				   ;  v d1 d2 d3 d4
					) collect
			(setf v (v+ v e)))))
	  (defparameter *vs* vs)
	  (store-vertex-list vs)
	  (progn
	    
	   #+nil (unless *initialized*
	     (vertex-attrib-pointer 0 3 :float nil 0 *vertex-buf*)
	     (enable-vertex-attrib-array 0)
	     (setf *initialized* t))
	   
	   (#_glVertexAttribPointer 0 3 #$GL_FLOAT 0 0 *vertex-buf*)

		
	    (color .9 .3 .3 1)
	    ;	    (draw-arrays :lines 0 2)
	    ))))))
#+nil(defparameter *bdfa* (%int-to-ptr (%address-of *vertex-buf*)))
(paref *vertex-buf* :float 0)
#+nil
(defparameter egl::*draw-function* #'bla)

#+nil
(defparameter egl::*run-gl* nil)
#+nil
(defparameter egl::*run-gl* t)

#+nil
(defparameter egl::*draw-function*
  #'(lambda ()
      (clear-color 1 1 0 .1)
       (clear :color-buffer-bit)))

#+nil
(bordeaux-threads:make-thread #'egl:run-egl :name "gl")


#+nil
(let ((var 0)
      (var2 0))
  (defun draw (w)
    ;;    (enable :line-smooth)   (blend-func :src-alpha :one-minus-src-alpha)
    (clear :color-buffer :depth-buffer)
    (enable :depth-test)
    (color 1 1 1)
    (progn
      (incf var 4)
      (when (< 360 var)
	(setf var 0)))
    (progn
      (incf var2 1.2)
      (when (< 360 var2)
	(setf var2 0)))
    (with-pushed-matrix
      (progn
	(unless *geometry-parameters*
	  (setf *geometry-parameters* (get-parameters-from-maxima)))
	(unless (find-if-not #'numberp *geometry-parameters*)
	  (destructuring-bind (ll gg d h1 a h2) *geometry-parameters*
	    (let* ((alpha (* 30 (sin (* pi (/ 180) var2))))
		(m-sys (rotation-matrix alpha (v  0 0 1)))
		(d1 (m* m-sys (make-vec gg)))
		(d2 (m* m-sys (make-vec (- d) (- h1))))
		(d3 (m* m-sys (make-vec (- (- gg (+ d a))) (+ h1 h2))))
		(d4 (m* m-sys (make-vec (- a) (- h2)))))
	  (with-pushed-matrix
	    (color 1 1 1)
	    (translate .3 .1 .7)
	    (rotate -90 1 0 0)
	    (let ((s .05))
	      (scale s s s))

	    (let ()
	    #+nil (with-primitive :line-strip
	       (vertex 0 0)
	       (vertex-v d1)
	       (vertex-v (v+ d1 d2))
	       (vertex-v (v+ (v+ d1 d2) d3))
	       (vertex-v (v+ (v+ d1 d2) (v+ d3 d4)))))

	    (progn
	      (enable :light0)
	      (color 1 1 1)
	      (let* ((r 3d0)
		     (x (- r)) (y r)
		     (m (rotation-matrix -90 (v 0 1)))
		     (m2 (rotation-matrix (+ (* .5 alpha))
					  (v 0 0 1)))
		     (p (m* m2 (m* m (make-vec x x))))
		     (q (m* m2 (m* m (make-vec x y))))
		     (r (m* m2 (m* m (make-vec y y))))
		     (s (m* m2 (m* m (make-vec y x))))
		     (n (normalize (vx (v- q p) (v- q r)))))
		(enable :lighting)
	#+nil	(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
		(line-width 3)
	#+nil	(with-primitive :lines
		  (vertex 0 0)
		  (vertex-v (v* 50d0 n)))
		(line-width 1)))
	    (progn
	      (enable :light0)
	      (color 1 1 1)
	      (let* ((r 3d0)
		     (x (- r)) (y r)
		     (m (rotation-matrix -90 (v 0 1))) ;; rotate mirror into plane
		     (angle1 (* -180 (/ pi)
				(acos (v. (v* -1d0 (normalize d1)) 
					  (normalize d2)))))
		     (m2 (rotation-matrix (+ alpha 180 (* .5 angle1))
					  (v 0 0 1)))
		     (p (v+ d1 (m* m2 (m* m (make-vec x x)))))
		     (q (v+ d1 (m* m2 (m* m (make-vec x y)))))
		     (r (v+ d1 (m* m2 (m* m (make-vec y y)))))
		     (s (v+ d1 (m* m2 (m* m (make-vec y x)))))
		     (n (normalize (vx (v- q p) (v- q r)))))
		(defparameter *angle1* angle1)
		(enable :lighting)
#+nil		(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
#+nil		(with-primitive :lines
		  (vertex-v d1)
		  (vertex-v (v+ d1 (v* 10d0 n))))))


	    (progn
	      (enable :light0)
	      (color 1 1 1)
	      (let* ((r 3d0)
		     (x (- r)) (y r)
		     (m (rotation-matrix -90 (v 0 1))) ;; rotate mirror into plane
		     (angle2 (* -180 (/ pi)
				(acos (v. (v* -1d0 (normalize d2)) 
					  (normalize d3)))))
		     (m2 (rotation-matrix (+ alpha
					   *angle1*
					   (* .5 angle2))
					  (v 0 0 1)))
		     (p (v+ (v+ d2 d1) (m* m2 (m* m (make-vec x x)))))
		     (q (v+ (v+ d2 d1) (m* m2 (m* m (make-vec x y)))))
		     (r (v+ (v+ d2 d1) (m* m2 (m* m (make-vec y y)))))
		     (s (v+ (v+ d2 d1) (m* m2 (m* m (make-vec y x)))))
		     (n (normalize (vx (v- q p) (v- q r)))))
		(defparameter *angle2* angle2)
		(enable :lighting)
		#+nil
		(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
		#+nil
		(with-primitive :lines
		  (vertex-v (v+ d2 d1))
		  (vertex-v (v+ (v+ d2 d1) (v* 10d0 n))))))
	    
	    (progn
	      (enable :light0)
	      (color 1 1 1)
	      (let* ((r 3d0)
		     (x (- r)) (y r)
		     (m (rotation-matrix -90 (v 0 1))) ;; rotate mirror into plane
		     (angle3 (* 180 (/ pi)
				(acos (v. (v* -1d0 (normalize d3)) 
					  (normalize d4)))))
		     (m2 (rotation-matrix (+ alpha
					   *angle1*
					   *angle2*
					   180
					   (* .5 angle3)
					   )
					  (v 0 0 1)))
		     (center (v+ (v+ d3 d2) d1))
		     (p (v+ center (m* m2 (m* m (make-vec x x)))))
		     (q (v+ center (m* m2 (m* m (make-vec x y)))))
		     (r (v+ center (m* m2 (m* m (make-vec y y)))))
		     (s (v+ center (m* m2 (m* m (make-vec y x)))))
		     (n (normalize (vx (v- q p) (v- q r)))))
		(enable :lighting)
		#+nil (with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
		#+nil (with-primitive :lines
		  (vertex-v center)
		  (vertex-v (v+ center (v* 10d0 n))))))))))))))



