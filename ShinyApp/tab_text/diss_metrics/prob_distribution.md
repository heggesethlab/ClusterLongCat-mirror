One method is to assume that each set of data can be precisely modeled by a probability distribution. The first approach used to measure the dissimilarity between sequences focuses on the state distribution within each sequence. An example of distances between probability distributions we implement in the Shiny App is the squared Chi-square distance. Let $p_{j|x}$ be the proportion of time spent in state $j$ in sequence $x$, and $p_j$ the overall proportion of time spent in state $j$, the squared **Chi-squared** distance is:

$$d_{\text{chi}}^{2}(x,y) = \sum_{j=1}^{J} \frac{(p_{j|x}-p_{j|y})^2}{p_j}, \;\;\; \text{where J = total number of states}$$

The original Chi-squared measure is sensitive to time spent in states, yet is insensitive to the order and timing of the states. In order to overcome this limitation by stressing the local transitions, an extension of Chi-squared distance that focuses on the state distribution within successive - possibly overlapping - periods of the sequence is also implemented. Let $p_{j|x_k}$ be the proportion of time spent in state $j$ in the subsequence of $x$ over period $k$, $p_{j|k}$ be the overall proportion of time in state $j$ in the $k^{th}$ interval, the period-dependent Chi-squared distance is: 

$$d_{\text{chi}}^{2}(x,y) = \sum_{k=1}^{K} \sum_{j=1}^{|\Sigma|} \frac{(p_{j|x_k}-p_{j|y_k})^2}{p_{j|k}}, \;\;\; \text{where K = total number of intervals}$$

Here, bandwidth defines the length of the intervals and overlap denotes whether intervals should overlap or not. Note that bandwidth must be even when overlap is true. The default step is set to the length of the longest sequence and the default overlap is set to false.
