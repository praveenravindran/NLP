# NLP

Built a trigram hidden Markov model to identify gene names in biological text where the task was to implement a probabilistic model and a decoder for finding the most likely tag sequence for new sentences.

ACHIEVED THE EXPECTED F-1 score of 42 for this assignment.

A labeled training data set gene.train is used to train the model , a labeled and unlabeled version of the development set, gene.key and gene.dev is used for development and gene.test is used for testing purposes.

"Pre-Viterbi" steps - 

/* Use the training file to generate counts for every word and classify certain words with certain classes into RARE, HASNUMBER and ALLCAPS and LASTCAP is done.
* Reads development file and based on the specific word, the word is associated with specific class we can calculate the emission parameter
*/

Viterbi algorithm is used for building the model and the steps are given below
/* Viterbi Algorithm  
 * 
 * Definitions:-
 * v - Tag sequence at current position
 * u - Tag sequence at current - 1 position
 * w - Tag sequence at current - 2 position
 * n - length of sentence
 * k - 1 through n 
 * xk - word at position k
 * piFunction pi (k,u,v) = max value of ( pi(k-1,w,u) * q(v|w,u) * e(xk | v) ) for all w's
 * bpFunction bp(k,u,v) = the tag at w which gives you the maxpi(k,u,v)
 * Steps:-
 * Define possible tag sequences at position -1 and 0 as * and from 1 through n all the other possible tag sequences( n is length of the sentence) 
 * Iterate from k= 1 through n and for all possible combinations of u and v 
 	* find pi(k,u,v)
 	* Set bp(k,u,v)
 * Set y(n-1) and y(n) = arg max(pi(n,u,v) * q(STOP|u,v)) 
 * Iterate  through k = n-2 to 1 and set the tag at position k as yk = bp(k+2,y(k+1),y(k+2))
 * Set this Tag Sequence as the best tag sequence for the sentence
 * */
# NLP
