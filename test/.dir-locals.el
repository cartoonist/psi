;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((c++-mode . ((eval . (progn
                        (dap-register-debug-template
                         "LLDB::Run::psi-tests::with-args"
                         (list :type "lldb-mi"
                               :request "launch"
                               :name "LLDB::Run::psi-tests::with-args"
                               :target "~/workspace/github.com/cartoonist/psi/build/Debug/test/psi-tests"
                               :arguments (list "[index][iterator]")
                               :showDevDebugOutput nil
                               :cwd nil))

                        (dap-register-debug-template
                         "LLDB::Run::psi-tests"
                         (list :type "lldb-mi"
                               :request "launch"
                               :name "LLDB::Run::psi-tests"
                               :target "~/workspace/github.com/cartoonist/psi/build/Debug/test/psi-tests"
                               :showDevDebugOutput nil
                               :cwd nil))
                        )))))
