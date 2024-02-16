          seed =  -1

       seqfile = testbed/small/common-data/frogs.txt
      Imapfile = testbed/small/common-data/frogs.Imap.txt
       outfile = testbed/small/38/out/out.txt
      mcmcfile = testbed/small/38/out/mcmc.txt

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1  0 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

    speciesmodelprior = 0  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  KENTUCKY  CUMBERLAND  NEW  KANAWHA
                    5 6 4 5
                   ((KENTUCKY, CUMBERLAND),(NEW, KANWAHA));
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 330  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.004   # invgamma(a, b) for theta
      tauprior = 3 0.002    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
    locusrate = 1 0 0 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 400
      sampfreq = 2
       nsample = 1500