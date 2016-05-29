package geneNamesInBiologicalText;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class EparameterCalculator {
	private HashMap<String,Integer>wordCount;
	private HashMap<String,Integer>oMap;
	private HashMap<String,Integer>iGeneMap;
	private Integer oSequenceCount;
	private Integer iGeneSequenceCount;
	private HashMap<String,Integer>bigramCount;
	private HashMap<String,Integer>triGramCount;
	private ArrayList<String>sentences;
	
	// Given a tag y and a word x the formulae used for ePArameter calculation is e(x/y) = Count(y -> x) / Count(y) 
	public Double calculateEparameter(String word, String sequence){
		if(word == null || sequence == null)
			return null;
		Integer yGivenXcount = 0;
		boolean hasNumber = false;
		boolean allCaps = true;
		boolean lastCap = false;
		for (int i = 0; i < word.length(); i++){
		    char c = word.charAt(i);        
		    int asciiVal= (int) c;
		    if(asciiVal > 47 && asciiVal < 58){
			    	hasNumber = true;
			} else if(asciiVal < 65 || asciiVal > 90){
		    	 allCaps = false;
			} else if(i == word.length() - 1){
		    	char last = word.charAt(i);
		    	int asciiLastVal= (int) last;
		    	if(asciiLastVal >= 65 && asciiLastVal <= 90 && allCaps == false){
		    		lastCap = true;
				}	
		    }
		}
		if(sequence.equals("o")){
			if(oMap.get(word) == null && iGeneMap.get(word) == null){
				if(hasNumber){
					yGivenXcount = oMap.get("_HASNUMBER_");
				}else if(allCaps){
					yGivenXcount = oMap.get("_ALLCAPS_");
				}else if(lastCap){
					yGivenXcount = oMap.get("_LASTCAP_"); 
				} else {
					yGivenXcount = oMap.get("_RARE_");
				}	
			} else if (oMap.get(word) != null ){
				yGivenXcount = oMap.get(word);
			} else{
				return 0.0;
			}
			double val = ((double)yGivenXcount / (double)this.oSequenceCount);
			return val;
		}
		if(sequence.equals("i-gene")){
			if(iGeneMap.get(word) == null && oMap.get(word) == null){
				if(hasNumber){
					yGivenXcount = this.iGeneMap.get("_HASNUMBER_");
				}else if(allCaps){
					yGivenXcount = this.iGeneMap.get("_ALLCAPS_");
				} else if(lastCap ){
					yGivenXcount = this.iGeneMap.get("_LASTCAP_");
				//	System.out.println("line 85 " + word + "---yGivenXcount---" + yGivenXcount );
				} else {
					yGivenXcount = this.iGeneMap.get("_RARE_");
				}
			} else if(iGeneMap.get(word) != null){
				yGivenXcount = this.iGeneMap.get(word);
			} else{
				return 0.0;
			}
			double val = ((double)yGivenXcount/(double)this.iGeneSequenceCount);
			return val;
		}
		else{
			return null;
		}
	}
	
	// This is a unigram tagger and works purely based on the tag of every word without considering what the tags of the previous words are.
	public void tagger(File file){
		FileReader fileReader;
		String line = null;
		String seq = null;
		List<String> text = new ArrayList<String>();
		
		CorpusProcessor corpusProcessor = new CorpusProcessor();
		this.wordCount=corpusProcessor.getWordCount();
		this.iGeneMap=corpusProcessor.getiGeneMap();
		this.oMap=corpusProcessor.getoMap();
		this.iGeneSequenceCount=corpusProcessor.getiGeneSequenceCount();
		this.oSequenceCount=corpusProcessor.getOSequenceCount();
		this.bigramCount=corpusProcessor.getBigramCount();
		this.triGramCount=corpusProcessor.getTriGramCount();
		this.sentences = new ArrayList<String>();
		
		StringBuilder sb = new StringBuilder();
		
		try {
			fileReader = new FileReader(file);
			BufferedReader bufferReader = new BufferedReader(fileReader);
			try {
				while((line = bufferReader.readLine()) != null){
					String testLine =line;
					String newLine = null;
					double oCheck=0.0;
					double iGeneCheck=0.0;
					if(line.trim().length() > 0 ){
						sb.append(line + " ");
						oCheck = calculateEparameter(testLine,"o");
						iGeneCheck = calculateEparameter(testLine,"i-gene");
						if(oCheck > iGeneCheck){
							seq ="O";
							newLine = line +" "+seq;
							text.add(newLine);
						} else if(iGeneCheck > oCheck){
							seq ="I-GENE";
							newLine = line +" "+seq;
							text.add(newLine);
						} else if(iGeneCheck == oCheck){
							seq ="I-GENE";
							newLine = line +" "+seq;
							text.add(newLine);
						}
					} else{
						this.sentences.add(sb.toString());
						sb = new StringBuilder();
						text.add(line);
					}
				}
				bufferReader.close();
				fileReader.close();
				
				Path devFile = null;			
				if(file.getName().equals("gene.dev")){
					 devFile = Paths.get("gene_dev.p1.out");
				} else {
					 devFile = Paths.get("gene_test.p1.out");
				}
				Files.write(devFile, text, Charset.forName("UTF-8"));
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public HashMap<String, Integer> getWordCount() {
		return wordCount;
	}

	public void setWordCount(HashMap<String, Integer> wordCount) {
		this.wordCount = wordCount;
	}

	public HashMap<String, Integer> getoMap() {
		return oMap;
	}

	public void setoMap(HashMap<String, Integer> oMap) {
		this.oMap = oMap;
	}

	public HashMap<String, Integer> getiGeneMap() {
		return iGeneMap;
	}

	public void setiGeneMap(HashMap<String, Integer> iGeneMap) {
		this.iGeneMap = iGeneMap;
	}

	public Integer getoSequenceCount() {
		return oSequenceCount;
	}

	public void setoSequenceCount(Integer oSequenceCount) {
		this.oSequenceCount = oSequenceCount;
	}

	public ArrayList<String> getSentences() {
		return sentences;
	}

	public void setSentences(ArrayList<String> sentences) {
		this.sentences = sentences;
	}

	public Integer getiGeneSequenceCount() {
		return iGeneSequenceCount;
	}

	public void setiGeneSequenceCount(Integer iGeneSequenceCount) {
		this.iGeneSequenceCount = iGeneSequenceCount;
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
