          seed =  -1

       seqfile = Users/emily/Desktop/BPP/full_phylip.txt
      Imapfile = Users/emily/Desktop/BPP/kentucki.Imap.txt
       outfile = kentucki_A00.out.txt
      mcmcfile = kentucki_A00.mcmc.txt

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2 * speciesdelimitation algorithm1 finetune (a m)
*  speciesdelimitation = 1 1 2 0.5  * speciesdelimitation algorithm1 finetune (a m)
*          speciestree = 1 * species-tree by NNI

 *speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
 *speciesmodelprior has no effect for A00

  species&tree = 5  KENTUCKY  CUMBERLAND  NEW  KANAWHA	PETRAEUS
                    5 6 4 5
                   (((KENTUCKY, CUMBERLAND),(NEW, KANWAHA)), PETRAEUS);


       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 330    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.004 e   # gamma(a, b) for theta
      tauprior = 3 0.004   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

      finetune = 1: .01 .0001 .005 .0005 .2 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 8000
      sampfreq = 2 
       nsample = 100000
       threads = 2 1 1
