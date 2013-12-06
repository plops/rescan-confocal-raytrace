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


(eval-when (:execute :compile-toplevel :load-toplevel) (progn
							 (ql:quickload :cl-glut)
							 (ql:quickload :cl-opengl)
							 (ql:quickload :cl-glu)))

(defpackage :g (:use :cl :gl))
(in-package :g)


(deftype vec ()
  `(simple-array double-float (3)))

(defun make-vec (&optional (x 0d0) (y 0d0) (z 0d0))
  (declare ;(type double-float x y z)
	   (values vec &optional))
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
	   (values double-float &optional))
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
	      (values vec &optional))
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
	   (values vec &optional))
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
	   (values vec &optional))
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
	   (values double-float &optional))
  (let ((l2 (v. v v)))
    (declare (type double-float l2))
    (sqrt l2)))

(defun normalize (v)
  (declare (type vec v)
	   (values vec &optional))
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
	   (values mat &optional))
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
	   (values vec &optional))
  (let ((res (v)))
    (with-arrays (res v m)
      (dotimes (i 3)
	(dotimes (j 3)
	  (incf (res i) (* (m i j) (v j))))))
    res))
#+nil
(m* (rotation-matrix 90 (v 0d0 0d0 1d0)) (v 1d0))


(defclass planet-window (glut:window)
  ((year :accessor year :initform 0)
   (day :accessor day :initform 0))
  (:default-initargs
   :pos-x (- 1366 500) :pos-y 100 :width 500 :height 500
   :mode '(:double :rgb :depth) :title "planet.lisp"))

(defmethod glut:display-window :before ((w planet-window))
  (clear-color 0 0 0 0)
  (shade-model :flat))



(defun rotate-around-axis (x y z theta rx ry rz)
  (declare (type double-float x y z theta rx ry rz)
	   (values double-float double-float double-float &optional)
	   (optimize (speed 3) (debug 1) (safety 1)))
  (let* ((c (cos theta))
	 (s (sin theta))
	 (_c (- 1 c)))
    (values (+ (* (+ c (* _c rx rx)) x)
	       (* (- (* _c rx ry) (* s rz)) y)
	       (* (+ (* _c rx rz) (* s ry)) z))
	    (+ (* (+ (* _c rx ry) (* s rz)) x)
	       (* (+ c (* _c ry ry)) y)
	       (* (- (* _c ry rz) (* s ry)) z))
	    (+ (* (- (* _c rx rz) (* s ry)) x)
	       (* (+ (* _c ry rz) (* s rx)) y)
	       (* (+ c (* _c rz rz)) z)))))

(defmacro rvertex (x y z theta rx ry rz)
  `(let ((l (expt (+ (* rx rx)
		     (* ry ry)
		     (* rz rz)) -.5d0))) 
     (multiple-value-bind (a b c) (rotate-around-axis (* 1d0 ,x) (* 1d0 ,y) (* 1d0 ,z) 
						      (* pi (/ 180d0) ,theta) 
						      (* l ,rx) (* l ,ry) (* l ,rz))
       (vertex a b c))))

(defun reflect (x y z nx ny nz)
  (declare (type double-float x y z nx ny nz)
	   (values double-float double-float double-float &optional)
	   (optimize (speed 3) (debug 1) (safety 1)))
  (let ((-2kn (* -2 (+ (* x nx) (* y ny) (* z nz)))))
    (values (+ x (* nx -2kn))
	    (+ y (* ny -2kn))
	    (+ z (* nz -2kn)))))

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
   (sb-ext:run-program "/usr/local/bin/maxima" 
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

(defvar *geometry-parameters* nil)

(defvar *angle1* 0)
(defvar *angle2* 0)

(let ((var 0)
      (var2 0))
  (defun draw (w)
    (enable :line-smooth)
    (blend-func :src-alpha :one-minus-src-alpha)
					;   (reshape w (planet-window-wi))
    (clear-color .4 .5 .3 0)
    (matrix-mode :modelview)
    (load-identity)
    (glu:look-at 7 (+ 5 (* .2 (+ (* .9 (sin (* pi (/ 180) var2))) (sin (* pi (/ 180) var))))) 2 (* .1 (sin (* pi (/ 180) var2)) #+nil (expt (random 1.0) 4)) 0 (* 0 (random .02)) 0 1 0)
    #+nil (glu:look-at 1 8 1  
		 0 0 0
		 0 1 0)
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
     (rotate (+ 20 (* 1.8 (+ (* .1 (sin (* pi (/ 180d0) var))) (sin (* pi (/ 180d0) var2))))) 0 1 0)
     (line-width 3)
     (with-primitive :lines
       (color 1 0 0) (vertex 0 0) (vertex 2 0)
       (color 0 1 0) (vertex 0 0) (vertex 0 2)
       (color 0 0 1) (vertex 0 0) (vertex 0 0 2))
     (line-width 1)
     (with-primitive :lines
       (color .6 0 0) (vertex 0 0) (vertex 10 0)
       (color 0 .6 0) (vertex 0 0) (vertex 0 10)
       (color 0 0 .6) (vertex 0 0) (vertex 0 0 10))
     (color .2 .2 .1)
     (with-primitive :lines
       (loop for i from 1 below 10 do (vertex 0 i) (vertex 10 i))
       (loop for i from 1 below 10 do (vertex 0 0 i) (vertex 0 10 i))
       (loop for i from 1 below 10 do (vertex i 0 0) (vertex i 0 10)))
     (color 1 1 1)
    
     ;(rotate 90 0 0 1)
     
     
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

	    (with-primitive :lines
	      (color 1 0 0) (vertex 0 0) (vertex 20 0)
	      (color 0 1 0) (vertex 0 0) (vertex 0 20)
	      (color 0 0 1) (vertex 0 0) (vertex 0 0 20))
	    (color 1 1 1)

	    ;(enable :line-stipple)
	    (line-stipple 2 #b0001110001000111)
	    (let ()
	     (with-primitive :line-strip
	       (vertex 0 0)
	       (vertex-v d1)
	       (vertex-v (v+ d1 d2))
	       (vertex-v (v+ (v+ d1 d2) d3))
	       (vertex-v (v+ (v+ d1 d2) (v+ d3 d4)))))
	    (disable :line-stipple)
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
		(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
		(line-width 3)
		(with-primitive :lines
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
		(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
		(with-primitive :lines
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
		(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
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
		(with-primitive :quads
		  (normal-v n)
		  (loop for e in (list p q r s) do
		       (vertex-v e)))
		(disable :lighting)
		(color 1 1 0)
		(with-primitive :lines
		  (vertex-v center)
		  (vertex-v (v+ center (v* 10d0 n)))))))))))
     
     (rotate (+ var (year w)) 0 1 0)
     (translate 3 0 0)
     (rotate (+ (* 8 var) (day w)) 0 1 0))
   (glut:swap-buffers)))




(defun vertex-v (v)
  (declare (type vec v))
  (with-arrays (v)
    (vertex (v 0) (v 1) (v 2)))
  v)

(defun normal-v (v)
  (declare (type vec v))
  (with-arrays (v)
    (normal (v 0) (v 1) (v 2)))
  v)




(defmethod glut:display ((w planet-window))
  (draw w)
  (sleep (/ 62))
  (glut:post-redisplay))

(defmethod glut:reshape ((w planet-window) width height)
  (viewport 0 0 width height)
  (matrix-mode :projection)
  (load-identity)
  (glu:perspective 60 (/ width height) 1 20)
  (matrix-mode :modelview)
  (load-identity)
  (glu:look-at 7 5 2 0 0 0 0 1 0))

(defmethod glut:keyboard ((w planet-window) key x y)
  (declare (ignore x y))
  (flet ((update (slot n)
           (setf (slot-value w slot) (mod (+ (slot-value w slot) n) 360))
           (glut:post-redisplay)))
    (case key
      (#\d (update 'day 10))
      (#\D (update 'day -10))
      (#\y (update 'year 5))
      (#\Y (update 'year -5))
      (#\Esc (glut:destroy-current-window)))))

(defun rb-glut ()
  (glut:display-window (make-instance 'planet-window)
   ))


#+nil
(sb-thread:make-thread #'(lambda () (rb-glut))
		       :name "gl")


