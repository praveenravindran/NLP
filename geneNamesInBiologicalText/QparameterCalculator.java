package geneNamesInBiologicalText;


import java.util.HashMap;


public class QparameterCalculator {
	private HashMap<String,Integer>bigramCount;
	private HashMap<String,Integer>triGramCount;
	
	/*Given the current word tag yi and the previous word tags yi-2 and y i-1 the q parameter is 
	q(yi|y(i-2),y(i-1)) = count(y(i-2),y(i-1),y(i)) / count (y(i-2),y(i-1))*/
	public double calculateQparameter(String yiMinusTwo,String yiMinusOne, String yi){
		String triGram = yiMinusTwo+","+yiMinusOne+","+yi;
		String biGram = yiMinusTwo+","+yiMinusOne;
		Integer triGramValue = this.triGramCount.get(triGram);
		Integer biGramValue = this.bigramCount.get(biGram);
		double qParameter = (double)triGramValue/(double)biGramValue;
		return qParameter;
	}
	
	public HashMap<String, Integer> getBigramCount() {
		return bigramCount;
	}

	public void setBigramCount(HashMap<String, Integer> bigramCount) {
		this.bigramCount = bigramCount;
	}

	public HashMap<String, Integer> getTriGramCount() {
		return triGramCount;
	}

	public void setTriGramCount(HashMap<String, Integer> triGramCount) {
		this.triGramCount = triGramCount;
	}
}
