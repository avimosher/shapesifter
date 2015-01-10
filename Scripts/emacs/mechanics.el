;#####################################################################
; PhysBAM Syntax Highlighting
;#####################################################################

(setq font-lock-keyword-case-fold-search nil) ; Need to be case sensitive
(font-lock-add-keywords 'c++-mode  '(("[^[:lower:]]\\([[:upper:]][[:upper:][:digit:]_]*\\)[<> ,]" 1 font-lock-type-face t))) ; match class names

(make-face 'font-lock-operators-face)
(font-lock-add-keywords 'c++-mode '(("[;<>:=!]\\|->" . 'font-lock-operators-face)))

(make-face 'font-lock-preprocessor-face)
(font-lock-add-keywords 'c++-mode '(("\\#[a-zA-Z0-9]*[ ]" . 'font-lock-preprocessor-face)))


;#####################################################################
; PhysBAM Helper Routines
;#####################################################################

(defun physbam-reduce (f x)
  (if (eq (cdr x) nil)
      (car x)
    (funcall f (car x) (physbam-reduce f (cdr x))))) 

(defun physbam-filter (f list)
  (let (filtered)
    (dolist (x list filtered)
      (if (funcall f x) (setq filtered (cons x filtered)) nil))))

;#####################################################################
; PhysBAM Navigation Commands
;#####################################################################

(defun physbam-header-flip ()
  "Find the .h file for this .C file (or vice versa)."
  (interactive)
  (let ((dotc (string-match "[.]\\(cpp\\|c\\)$" (buffer-file-name)))
        (doth (string-match "[.]h$" (buffer-file-name))))
    (if dotc
        (find-file (concat (substring (buffer-file-name) 0 dotc) ".h"))
      (if doth
          (find-file (concat (substring (buffer-file-name) 0 doth) ".cpp"))
        (message "Not a cpp or h file!!")))))

(defun physbam-dimension-flip ()
  "Find the .h file for this .C file (or vice versa)."
  (interactive)
  (let ((buf (buffer-name)))
    (if (string-match "\\(.+\\)\\(1D\\|2D\\|3D\\)\\(.+\\)"  buf)
      (let ((current (match-string 2 buf))
            (left (match-string 1 buf))
            (right (match-string 3 buf)))
        (find-file (concat left
                (cond ((string= "1D" current) "2D") ((string= "2D" current) "3D") ((string= "3D" current) "1D"))
                right))))))

(defun physbam-open-parent ()
  "Open header of parent class"
  (interactive)
  (let ((doth (string-match "[.]h$" (buffer-file-name))))
    (if doth
        (let ((save_point (point)))
          (goto-char 0)
          (if (re-search-forward ":\\(?:public\\|protected\\|private\\) \\(\\(?:\\w\\|_\\)*\\)" nil t)
              (let ((filename (concat (match-string 1) ".h")))
                (goto-char save_point)
                (if (file-exists-p filename)
                    (find-file filename)
                  (let ((physbam_filename (find-physbam-file filename)))
                    (if (file-exists-p physbam_filename)
                        (find-file physbam_filename)
                      (message "Can't find header of parent")))))
            (goto-char save_point)
            (message "No parent class")))
      (message "Not a header file"))))

(defun mechanics-grep-cpp-and-headers (directory_prefix querystr)
  (let ((path (concat (getenv "MECHANICS") "/" directory_prefix "/")))
  (grep (concat "(cd " path ";find -name '*.h' -o -name '*.cpp' | xargs grep -n -e \"" querystr "\") | sed 's@^@" path "@' # "))))

(defun mechanics-grep-library (querystr)
  "grep in Library"
  (interactive "sGrep $MECHANICS/Library for:")
  (mechanics-grep-cpp-and-headers "Library" querystr))


(defun mechanics-fix-includes ()
  "Fix includes in current file"
  (interactive)
  (save-buffer)
  (shell-command (format "%s/Scripts/misc/fix_headers_file.sh %s" (getenv "MECHANICS") (buffer-file-name)))
  (print "Test2")
  (revert-buffer-no-prompt))

(defun physbam-insert-header ()
  (interactive)
  (physbam-insert-copyright)
  (let ((classname (substring (buffer-name) 0 -2)))
    (insert (format "#ifndef __%s__\n" classname))
    (insert (format "#define __%s__\n\n" classname))
    (insert "namespace PhysBAM{\n\n")
    (insert "template<class T>\n")
    (insert (format "class %s\n" classname))
    (insert "{\n")
    (insert "public:\n\n")
    (insert "//#####################################################################\n")
    (insert "};\n")
    (insert "}\n")
    (insert "#endif\n")))

(defun physbam-make-cpp-from-header ()
  (interactive)
  (goto-char 0)
  (goto-line 7)
  (copy-region-as-kill (point-min) (point))
  (physbam-header-flip)
  (yank)
  (let ((classname (substring (buffer-name) 0 -4)))
    (insert (format "#include \"%s.h\"\n" classname))
    (insert "using namespace PhysBAM;\n")
    (insert (format "template class %s<float>;\n" classname))
    (insert (format "template class %s<double>;\n" classname))))

;#####################################################################
; PhysBAM Build / Run Commands
;#####################################################################

(defun mechanics-setup-compile-command (write_settings)
  (setq compile-command (format "scons --warn=no-deprecated -Q --implicit-cache -D TYPE=debug -j 4 CWD=%s" (getenv "MECHANICS")))
  (message (format "New compile command is: %s" compile-command))
  (setq compilation-read-command nil))

(defun mechanics-compile ()
  (interactive)
  (call-interactively 'compile))


;(setq physbam-project-type "release")
;(setq physbam-compile-count 4)

;(setq compile-command (format "make -k" physbam-project-directory))

;(setq tags-file-name (format "%s/TAGS" (getenv "PHYSBAM")))
;(setq physbam-compile-mode "single")
; Setup default project to be currenct directory
; Read project settings
;(physbam-read-project-settings)
; NOTE ALL FUNCTIONS THAT MODIFY STATUS AND SAVE SHOULD BE BELOW ABOVE READ PROJECT SETTINGS
; (physbam-set-compiler (if (or (string= (getenv "PLATFORM") "opteron") (string= (getenv "PLATFORM") "nocona"))  "gcc-4.0.1-64" "gcc-4.0.1"))
;(mechanics-set-compiler "g++")
(mechanics-setup-compile-command nil)
(setq truncate-partial-width-windows nil)
;(setq compilation-scroll-output t) ; scroll to end by default

(defconst physbam-c-style
  '((c-basic-offset . 4)
        (c-offsets-alist . ((string . -1000)
                        (c . c-lineup-C-comments)
                        (defun-open . 0)
                        (defun-close . 0)
                        (defun-block-intro . +)
                        (class-open . 0)
                        (class-close . 0)
                        (inline-open . 0)
                        (inline-close . 0)
                        (func-decl-cont . +)
                        (knr-argdecl-intro . 5)
                        (knr-argdecl . 0)
                        (topmost-intro . 0)
                        (topmost-intro-cont . 0)
                        (member-init-intro . +)
                        (member-init-cont . -1)
                        (inher-intro . +)
                        (inher-cont . c-lineup-multi-inher)
                        (block-open . 0)
                        (block-close . 0)
                        (brace-list-open . 0)
                        (brace-list-close . 0)
                        (brace-list-intro . +)
                        (brace-list-entry . 0)
                        (statement . 0)
                        (statement-cont . +)
                        (statement-block-intro . +)
                        (statement-case-intro . +)
                        (statement-case-open . +)
                        (substatement . +)
                        (substatement-open . 0)
                        (case-label . +)
                        (access-label . -)
                        (label . *)
                        (do-while-closure . 0)
                        (else-clause . 0)
                        (comment-intro . c-lineup-comment)
                        (arglist-intro . +)
                        (arglist-cont . 0)
                        (arglist-cont-nonempty . +)
                        (arglist-close . 0)
                        (stream-op . c-lineup-streamop)
                        (inclass . +)
                        (cpp-macro . -1000)
                        (friend . 0)
                        (extern-lang-open . 0)
                        (extern-lang-close . 0)
                        (inextern-lang . +)
			(namespace-open . 0)
			(namespace-close . 0)
			(innamespace . 0)
                        (template-args-cont . +)))    
    (c-hanging-braces-alist . ((class-open before after)
                               (class-close before)
                               (defun-open before after)
                               (defun-close before after)
                               (inline-open before)
                               (inline-close after)
                               (brace-list-open)
                               (brace-list-close)
                               (brace-list-intro)
                               (brace-list-entry)
                               (block-open after)
                               (block-close)
                               (substatement-open after)
                               (substatement-case-open)
			       (namespace-open after)
			       (namespace-close before after)
                               (extern-lang-open after)
                               (extern-lang-close after)))
    (c-hanging-colons-alist . ((case-label)
                               (label)
                               (access-label after)
                               (member-init-intro before)
                               (inher-intro)))
    (c-cleanup-list         . ((defun-close-semi)
                               (scope-operator)))
    (c-hanging-semi&comma-criteria . (c-semi&comma-no-newlines-before-nonblanks))
    (c-echo-syntactic-information-p . t))
  "PhysBAM C++ Style")

(c-add-style "physbam" physbam-c-style)
