;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((c++-mode . ((eval . (progn
                        (dap-register-debug-template
                         "LLDB::Run::rcrs_benchmark"
                         (list :type "lldb-mi"
                               :request "launch"
                               :name "LLDB::Run::rcrs_benchmark"
                               :target "~/workspace/github.com/cartoonist/psi/build/Debug/tools/rcrs_benchmark"
                               :arguments (list "--threads=1")
                               :showDevDebugOutput nil
                               :cwd nil))
