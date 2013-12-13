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


(eval-when (:execute :compile-toplevel :load-toplevel) 
  (progn
    (ql:quickload :cl-glut)
    (ql:quickload :external-program)
    (ql:quickload :cl-opengl)
    (ql:quickload :cl-glu)
    (ql:quickload :bordeaux-threads)
    (ql:quickload :alexandria)))


(defpackage :g (:use :cl :gl))
(in-package :g)


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


(defclass planet-window (glut:window)
  ((year :accessor year :initform 0)
   (day :accessor day :initform 0))
  (:default-initargs
   :pos-x (- 1366 500) :pos-y 100 :width 500 :height 500
   :mode '(:double :rgb :depth) :title "planet.lisp"))

(defmethod glut:display-window :before ((w planet-window))
  (clear-color 0 0 0 0)
  (shade-model :flat))

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
solution:[f1,f2,ll,gg,d,sol[1][3],a,sol[1][4]];
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


(defvar *geometry-parameters* (list 50 100
				    300
				    75.0
				    7.5
				    52.492676
				    15.0
				    52.492676))

(defvar *angle1* 0)
(defvar *angle2* 0)

(defun draw-mirror (center angle &key (normal-len 10d0))
  (declare (type vec center))
  (enable :light0)
  (color 1 1 1)
  (let* ((r 3d0)
	 (x (- r)) (y r)
	 (m (rotation-matrix -90 (v 0 1)))
	 (m2 (rotation-matrix angle
			      (v 0 0 1)))
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
    (line-width 3)
    (with-primitive :lines
      (vertex-v center)
      (vertex-v (v+ center (v* normal-len n))))
    (line-width 1))
  angle)

(defun draw-lens (center angle &key (normal-len 10d0))
  (declare (type vec center))
  (draw-mirror center angle :normal-len normal-len)
  angle)

(defclass disk ()
  ((center :accessor disk-center :initarg :center 
	   :initform (alexandria:required-argument)
	   :type vec)
   (normal :accessor disk-normal :initarg :normal
	   :initform (alexandria:required-argument)
	   :type vec)
   (radius :accessor disk-radius :initarg :radius
	  :initform (alexandria:required-argument)
	  :type double-float)))

(defclass lens (disk)
  ((focal-length :accessor focal-length :initarg :focal-length
		 :initform (alexandria:required-argument)
		 :type double-float)))

(defmethod print-object ((lens lens) stream)
  (with-slots (focal-length) lens
   (format stream "#<lens f: ~2,2f>" focal-length)))

(defclass ray ()
  ((start :accessor start :initarg :start
	  :initform (alexandria:required-argument)
	  :type vec)
   (direction :accessor direction :initarg :direction
	  :initform (alexandria:required-argument)
	  :type vec)))

(defmethod print-object ((ray ray) stream)
  (with-slots (start direction) ray
   (format stream "#<ray start: ~a dir: ~a>" start direction)))

(defgeneric intersect (ray object))

(defmethod intersect ((ray ray) (disk disk))
  "Find the point where a ray intersects a plane."
  #+sbcl (declare (values vec &optional))
  (with-slots (center normal) disk
    (with-slots (start direction) ray
      (let* ((hess-dist (v. center normal)) ;; distance of plane to origin
	     (eta (/ (- hess-dist (v. normal start)) ;; scaling along ray to hit exactly on plane
		   (v. normal direction))))
	(v+ start (v* eta direction))))))

#+nil
(intersect 
 (make-instance 'ray :start (v 0 1 -1) :direction (v 0 0 1))
 (make-instance 'plane :normal (v 0 0 1) :center (v)))


(define-condition ray-lost () ())

(defgeneric refract (ray object))

;;
;;		   --\
;;		  |   -----\   intersection
;;		  | 	    --+
;;	   |	   \	   ---|--\
;;	   |	   | -----/   |   ----\
;;	   |   ----+/  	r'    |	       ----\
;;       f |--/    c	      |		    ----\      dir
;;         |  	r	      |		         ---\
;;	   |  		      |rho		     ----\
;;	   |   		      |				  ----\
;;	   |   		      |				       ----\
;;	   |   		      |				  phi	    ----\
;;  -------+----------<-------+-------------------------------------------
;;	   |             n    |	center
;;	   |     	      |
;;	   |     	      |
;;	   |    	      |
;;	      		      |
;;	     	      f	      |

(defmethod refract ((ray ray) (lens lens))
  "Return new ray after refraction on thin lens. In general you will
have to normalize its direction. The refraction on an objective needs
the non-normalized result. When the ray doesn't hit the lens the
condition RAY-LOST is signalled."
  (declare (values ray &optional))
  (with-slots (start direction) ray
    (with-slots (center normal focal-length radius) lens
      (assert (< (abs (- 1 (norm direction))) 1e-12))
      (assert (< (abs (- 1 (norm normal))) 1e-12))
      (let* ((intersection (intersect ray lens))
	     (rho (v- intersection center))
	     (cosphi (v. normal direction)))
	(when (< radius (norm rho))
	  (error 'ray-lost))
	(make-instance 'ray
		       :start intersection
		       :direction (v- (v* (/ focal-length cosphi) direction)
				      rho))))))
#+nil
(handler-case 
    (refract (make-instance 'ray :start (v 0 .1 -10)
			    :direction (v 0 0 1))
	     (make-instance 'lens 
			    :focal-length 10.0
			    :center (v)
			    :normal (v 0 0 1)
			    :lens-radius .2))
  (ray-lost () nil))

;;   sketch of the mirror for incoming parallel light
;;   --------------------+-----------------------
;; 		      /|\
;; 		     / | \
;; 		    / n|  \	       N=n*(- (p.n))
;; 		q  /   |   \  p	       p+N=r
;; 	          /    v    \	       q=N+r
;; 	         /           \
;; 	        /      |      \	       q=p-2(p.n)*n
;; 	       /       |       \
;; 	      /       N|        \
;; 	     /         |         \
;; 	    /          |     r    \
;; 	   /   	       v<----------\
;; p .. ray-direction
;; N .. mirror-normal

(defgeneric reflect (ray object))

(defmethod reflect ((ray ray) (disk disk))
  "Return reflected ray. If the ray isn't inside of the radius return
signal RAY-LOST."
  (declare (values ray &optional))
  (with-slots (start direction) ray
   (with-slots (center normal radius) disk
     (assert (< (abs (- 1 (norm normal))) 1e-12))
     (assert (< (abs (- 1 (norm direction))) 1e-12))
     (let ((intersection (intersect ray disk)))
       (when (< (norm (v- intersection center)) radius)
	 (signal 'ray-lost))
       (let ((dir (v+ direction (v* (* -2d0 (v. direction normal))
					normal))))
	 (make-instance 'ray 
			:start intersection 
			:direction (normalize dir)))))))

#+nil
(reflect 
 (make-instance 'ray :start (v 0 0 10)
		:direction (v 0 0 -1))
 (make-instance 'disk :radius 1.0
		:center (v)
		:normal (v 0 0 1)))

(let ((var 0)
      (var2 0))
  (defun draw (w)
    (enable :line-smooth)
    (blend-func :src-alpha :one-minus-src-alpha)
					;   (reshape w (planet-window-wi))
    (clear-color .4 .5 .3 0)
    (matrix-mode :modelview)
    (load-identity)
    (glu:look-at 
     7 
     (+ 5 (* .2 (+ (* .9 (sin (* pi (/ 180) var2)))
		   (sin (* pi (/ 180) var))))) 
     2 
     
     (* .1 (sin (* pi (/ 180) var2)))            0        (* 0 (random .02)) 

     0 1 0)
    #+nil (glu:look-at 1 8 1  
		 0 0 0
		 0 1 0)
    (clear :color-buffer :depth-buffer)
   (enable :depth-test :normalize)
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
         
     (progn
      (unless *geometry-parameters*
	(setf *geometry-parameters* (get-parameters-from-maxima)))
      (unless (find-if-not #'numberp *geometry-parameters*)
       (destructuring-bind (f1 f2 ll gg d h1 a h2) *geometry-parameters*
	 (declare (ignorable ll))
	 (let* ((alpha (+ 10 (* 0 (sin (* pi (/ 180) var2)))))
		(m-sys (rotation-matrix alpha (v  0 0 1)))
		(d1 (m* m-sys (make-vec gg)))
		(d2 (m* m-sys (make-vec (- d) (- h1))))
		(d3 (m* m-sys (make-vec (- (- gg (+ d a))) (+ h1 h2))))
		(d4 (m* m-sys (make-vec (- a) (- h2)))))
	  (with-pushed-matrix
	    (translate .3 .1 .7)
	    (rotate -90 1 0 0)
	    (let ((s .05))
	      (scale s s s))

	    (with-primitive :lines ;; local coordinate system
	      (color 1 0 0) (vertex 0 0) (vertex 20 0)
	      (color 0 1 0) (vertex 0 0) (vertex 0 20)
	      (color 0 0 1) (vertex 0 0) (vertex 0 0 20))
	    
	    (with-primitive :line-strip  ;; optical axis
	      (color 1 1 1) (vertex 0 0) (vertex-v d1)
	      (vertex-v (v+ d1 d2)) (vertex-v (v+ (v+ d1 d2) d3))
	      (vertex-v (v+ (v+ d1 d2) (v+ d3 d4))))

	    (let ((angle1 (* -180 (/ pi)
			     (acos (v. (v* -1d0 (normalize d1)) 
				       (normalize d2)))))
		  (angle2 (* -180 (/ pi)
				(acos (v. (v* -1d0 (normalize d2)) 
					  (normalize d3)))))
		  (angle3 (* 180 (/ pi)
			     (acos (v. (v* -1d0 (normalize d3)) 
				       (normalize d4))))))
	      (draw-mirror (v) 
			   (+ (* .5 alpha) (* .5 angle1))
			   :normal-len 50d0)
	      (draw-mirror  d1
			    (+ alpha
			       180
			       (* .5 angle1)))
	      (draw-mirror  (v+ d2 d1) 
			    (+ alpha
			       angle1
			       (* .5 angle2)))
	      (draw-mirror  (v+ (v+ d3 d2) d1)
			    (+ alpha angle1 angle2 
			       180
			       (* .5 angle3)))
	      
	      (multiple-value-bind (center angle)
		  (get-point-along-polygon (list d1 d2 d3 d4) (* .99 ll .5 (+ 1 (sin (* var2 (/ 180) pi)))))
		(draw-lens center angle))
	      
	      (defparameter *bla* (list d1 d2 d3 d4))

#+nil	      (draw-lens)))))))
     
     (rotate (+ var (year w)) 0 1 0)
     (translate 3 0 0)
     (rotate (+ (* 8 var) (day w)) 0 1 0))
   (glut:swap-buffers)))

(loop for e in *bla* collect (* 180 (/ pi) (atan (aref e 1) (aref e 0))))

(defun get-point-along-polygon (rvs l &key (start (v)))
  "Given a polygon by relative vectors, find the point of the polygon
at a circumference length l. Return this point and the current
direction vector."
  (let ((len 0))
   (loop for p in rvs and index from 0 do
	(incf len (norm p))
	(when (< l len)
	  (return-from get-point-along-polygon
	    (values (v+ start (v* (- (norm p) (- len l)) (normalize p)))
		    (* -180 (/ pi) (atan (aref p 1) (aref p 0))))))
	(setf start (v+ start p)
	      p-old p))
   (error "requested length not within circumference of given polygon.")))

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
  (glu:perspective 60 (/ width height) 1 20))

(defmethod glut:keyboard ((w planet-window) key x y)
  (declare (ignore x y))
  (case key
    (#\Esc (glut:destroy-current-window))))

(defun rb-glut ()
  (glut:display-window (make-instance 'planet-window)))

#+nil
(bordeaux-threads:make-thread #'(lambda () (rb-glut)) :name "gl")



