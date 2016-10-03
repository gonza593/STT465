# STT465: Bayesian Statistical Methods (MSU)


* **Instructor:** Gustavo de los Campos ( gustavoc@msu.edu )
* **Time/Place:** MW 10:20am-11:40am A120 Wells Hall (WH)   
* **[Syllabus](https://github.com/gdlc/STT465/blob/master/STT465_Syllabus.pdf)**
* **[Required textbook](http://www.stat.washington.edu/people/pdhoff/book.php)**
* **[R-software](http://www.r-project.org/)**
* **[Tentative Schedule](https://github.com/gdlc/STT465/blob/master/SCHEDULE.pdf)**

## Homework
 * **[HW 1](https://www.dropbox.com/s/xj5uf9540udtdje/HW1_STT465.pdf?dl=0)** (Due W Sept 14 in class).
 * **HW 2**: Problems 3.3, 3.7 and 3.9 from the book (due  Wednesday, Sept 28th in class).
 
## Chapter 2. Belief, probability and exchangeability

 * Sections: 2.1-2.8, tpics covered:
  * Beliefs and properties of belief functions
  * Rules of probability
   * Total probability
   * Marginal probability
   * Bayes Rule
  * Probability functions as a way of describing beliefs
  * Independence
  * Conditional independence
  * Discrete RV (definition, pdf, probability of events, examples, includin Bernoulli, Binomial and Poisson)
  * Contionous RV (definition, CDF, pdf, examples: Normal)
  * Descriptors of Distributions: Expectation, Variance, quantiles.
  * Joint Marginal and Conditional Distributions (both for discrete and contionous)
  * The joint distribution of independent and IID RV
  * Exchangeability (definiton and example: Bernoulli)
  * Conditional independence (example: Bernoulli model)
  * De Finetty's Theorem (theorem and implications for building Bayesian models)
 
## Chapter 3. One Parameter Models
 * For the Beta-Binomial and Poisson-Gamma models we will cover:
    * Sampling Model
    * Maximum Likelihood Estimation
    * Large sample distribution of MLE
    * Large sample Frequentist CI (computation and interpretation)
    * Prior
    * Posterior
    * Bayesian inference
      * Maximum a posteriori
      * Posterior mean
      * Bayesian Credibility Regions

## Chapter 4. Monte Carlo Approximation
  * The method
  * Inference of arbitrary functions
  * The predictive distribution
    * In the Beta-Biomial model
    * In the Poisson-Gamma model

## Chapter 5. The Normal Model
  * The normal distribution: parameters, and properties
  * Likelihood for nomal model
    * MLE estimate of the mean and variance
    * Bias and variance of the MLE estimates
  * Bayesian model for the mean with known variance
    * Likelihood
    * Normal prior for the mean
    * Posterior distribution for the mean
    * The Bayesian estimate as a compromise between the information provided by the likelihood (MLE estiamate) and that conveyed by the prior.
    * The normal model for the mean and variance
      * Likelihood
      * The Scaled-invers Chi-square
      * Joint prior distribution for the mean and variance
      * Joint posterior
      * Alternative sampling schemes
## Chapter 6. The Gibbs Sampler
   * [Casella & George, AMSTAT, 1992](http://www.jstor.org/stable/2685208?Search=yes&resultItemClick=true&searchText=casella&searchText=the&searchText=gibbs&searchText=sampler&searchText=1992&searchUri=%2Faction%2FdoBasicSearch%3FQuery%3Dcasella%2Bthe%2Bgibbs%2Bsampler%2B1992%26amp%3Bacc%3Doff%26amp%3Bwc%3Don%26amp%3Bfc%3Doff%26amp%3Bgroup%3Dnone&seq=1#page_scan_tab_contents)
   * The Gibbs sample in the normal model for the mean and variance
