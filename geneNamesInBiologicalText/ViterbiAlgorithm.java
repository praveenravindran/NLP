package geneNamesInBiologicalText;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
 * Set this Tag Sequence
 * */

public class ViterbiAlgorithm {
	public ViterbiAlgorithm(){
		
	}
	private ArrayList<String>sentences;
	private HashMap<Integer,ArrayList<String>> possibleTagSeq;
	private HashMap<Integer,String> wordSequence;
	private Map<Integer,HashMap<String,Double>> piValues;
	private Map<Integer,HashMap<String,String>> backPointers;
	EparameterCalculator eparameterCalculator;
	QparameterCalculator qparameterCalculator;
	String sMinusOne = "*";
	String sMinusTwo = "*";
	List<String> lines = new ArrayList<String>();

	public double calculatePiVal(Integer k , String u , String v){
		double max = -10;
		double value = 0.0;
		String bp = "";
		String uv = u + v;
		if(k == 0 && u.equals("*") && v.equals("*"))
			return 1.0;
		else{
			ArrayList<String> wSequences = possibleTagSeq.get(k-2);
			String word = this.wordSequence.get(k);
			HashMap<String,Double> piMap = new HashMap<String,Double>();
			HashMap<String,String> bpMap = new HashMap<String,String>();
			if(this.piValues.get(k) != null){
				piMap = this.piValues.get(k);
			}
			if(this.backPointers.get(k) != null){
				bpMap = this.backPointers.get(k);
			}
			if(piMap.get(uv) == null){
				for(String w:wSequences){
					value = calculatePiVal(k-1,w,u) * this.qparameterCalculator.calculateQparameter(w, u, v) * this.eparameterCalculator.calculateEparameter(word, v);
					if(value > max){
						max = value;
						bp = w;
					}
				}
				bpMap.put(uv, bp);
				piMap.put(uv, max);
				this.piValues.put(k, piMap);
				this.backPointers.put(k, bpMap);
			} else {
				max = piMap.get(uv);
			}
		
		}
		return max;
	}
	
	public void doViterbi(String sentence){
		String words[] = sentence.split(" ");
		Integer sentenceLength = words.length - 1;
		ArrayList<String> uSequences = new ArrayList<String>();
		ArrayList<String> vSequences = new ArrayList<String>();
		String finalSequences[] = new String[sentenceLength + 1];
		this.possibleTagSeq = new HashMap<Integer,ArrayList<String>>();
		this.wordSequence = new HashMap<Integer,String>();
		ArrayList<String> sequences = new ArrayList<String>();
		double val = 0.0;
		sequences.add("*");
		possibleTagSeq.put(-1, sequences);
		possibleTagSeq.put(0, sequences);
		for(int j =1; j <= sentenceLength; j++){
			sequences = new ArrayList<String>();
			sequences.add("o");
			sequences.add("i-gene");
			possibleTagSeq.put(j, sequences);
			wordSequence.put(j, words[j]);
		}
		for(int k =1; k <= sentenceLength; k++){
			 uSequences = possibleTagSeq.get(k-1);
			 vSequences = possibleTagSeq.get(k);
			for(String u : uSequences){
				for(String v : vSequences){
					val = calculatePiVal(k,u,v);
				}
			}
		}
		//
		double maxProb = -10.0;
		finalSequences[0] = " ";
		for(String u:uSequences){
			for(String v:vSequences){
				double prob = calculatePiVal(sentenceLength,u,v) * this.qparameterCalculator.calculateQparameter(u, v, "stop");
				if(prob > maxProb){
					maxProb = prob;
					finalSequences[sentenceLength] = v;
					finalSequences[sentenceLength - 1] = u;
				}
			}
		}
		//
		for(int j = sentenceLength - 2; j >= 1; j-- ){
			HashMap<String,String> tempMap = this.backPointers.get(j+2);
			String bp = tempMap.get(finalSequences[j+1].toLowerCase()+finalSequences[j+2].toLowerCase());
			finalSequences[j] = bp;	
		}
		//
		for(int l =1 ; l <= sentenceLength; l++){
			lines.add(words[l] + " " + finalSequences[l].toUpperCase());
		}
	}
	
	public void helpViterbi(){
		File file2 = new File("gene.dev");
		this.eparameterCalculator = new EparameterCalculator();
		eparameterCalculator.tagger(file2);
		this.sentences = eparameterCalculator.getSentences();
		this.qparameterCalculator = new QparameterCalculator();
		qparameterCalculator.setBigramCount(eparameterCalculator.getBigramCount());
		qparameterCalculator.setTriGramCount(eparameterCalculator.getTriGramCount());
		String sentence = null;
		for (String sent:this.sentences) {
			sentence= " " + sent;
			piValues = new HashMap<Integer,HashMap<String,Double>>();
			backPointers = new HashMap<Integer,HashMap<String,String>>();
			this.doViterbi(sentence);
			this.lines.add(" ");
		}
		Path file = Paths.get("geneSeqOp.dev");
		Path oldDevFile = Paths.get("gene_dev.p1.out");
		try {
			Files.delete(oldDevFile);
			Files.write(file, lines, Charset.forName("UTF-8"));
			new File("geneSeqOp.dev").renameTo(new File("gene_dev.p1.out"));
			System.out.println("sccuess");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public static void main(String args[]){
		ViterbiAlgorithm viterbiAlgorithm = new ViterbiAlgorithm();
		viterbiAlgorithm.helpViterbi();
	}

}
