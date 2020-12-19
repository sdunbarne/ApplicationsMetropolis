# ApplicationsMetropolis
        For the Ising model the Metropolis algorithm takes a random walk
        through the configuration space, visiting more frequently
        occurring spin configurations more often.  Instead of choosing
        configurations randomly, then weighting them with \( \EulerE^{-E/kT}
        \), instead choose configurations with probability \( \EulerE^{-E/kT}
        \) and weight them evenly.
    
        For image reconstruction, the goal is to find the configuration
        maximizing \( \Prob{\omega \given \omega^{\text{blurred}}} \),
        called the maximum a posteriori estimate.  The technique is to
        formulate a new version of the Metropolis algorithm called \defn
        {Gibbs sampling}.
      
        Bayesian hierarchical models naturally describe the
        connections between data, observed parameters and other
        unobserved parameters, sometimes called \emph{latent
          variables}.  Suppose \( f(x_1, x_2, \dots, x_N) \) is a
        probability distribution in which the variables represent
        parameters of a statistical model.  The goal is to obtain
        point and interval estimates for these parameters. This fits
        into the Gibbs sampling framework, using Bayes' rule.
