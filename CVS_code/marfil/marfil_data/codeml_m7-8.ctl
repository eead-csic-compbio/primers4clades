      seqfile = seqfile.txt   * sequence data filename
     treefile = tree.txt      * tree structure file name
      outfile = results.txt   * main result file name

        noisy = 9      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = 0      * 0:user defined tree

      seqtype = 1      * 1:codons
    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61

        model = 0      * 0:one omega ratio for all branches

      NSsites = 7 8     * 0:one omega ratio (M0)
                       * 1:neutral (M1)
                       * 2:selection (M2)
                       * 3:discrete (M3)
                       * 7:beta (M7)
                       * 8:beta&w (M8)

        icode = 0      * 0:universal code

    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * initial or fixed kappa

    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated 
        omega = 5      * initial omega

                       *set ncatG for models M3, M7, and M8!!!
       *ncatG = 2      * # of site categories for M3 in Table 4
       ncatG = 10     * # of site categories for M7 and M8 in Table 4
