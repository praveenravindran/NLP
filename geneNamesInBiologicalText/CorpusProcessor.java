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
import java.util.Map;



public class CorpusProcessor {
	
	private HashMap<String,Integer> wordCount;
	private HashMap<String,Integer>oMap;
	private HashMap<String,Integer>iGeneMap;
	private HashMap<String,Integer>bigramCount;
	private HashMap<String,Integer>triGramCount;
	private Integer oSequenceCount;
	private Integer iGeneSequenceCount;
	
	public CorpusProcessor(){
		this.mapWords();
		this.replaceRareWords();
	}
	
	public void mapWords(){
		// Assign counts to all the words in gene.count
		String fileName = "gene.counts";
		String line = null;
		
		try{
			FileReader fileReader = new FileReader(fileName);
			BufferedReader bufferReader = new BufferedReader(fileReader);
			this.wordCount=new HashMap<String,Integer>();
			//this.rareWords=new HashMap<String,String>();
			this.oMap=new HashMap<String,Integer>();
			this.iGeneMap=new HashMap<String,Integer>();
			this.bigramCount=new HashMap<String,Integer>(); 
			this.triGramCount=new HashMap<String,Integer>();
			this.oSequenceCount=0;
			this.iGeneSequenceCount=0;
			
			while((line = bufferReader.readLine()) != null){
				String arr[] = line.split(" ");
				if(arr[1].toLowerCase().equals("wordtag")){
					if(arr[3] != null && arr[3].trim().length() > 0){
						String word = arr[3];
						if(this.wordCount.get(word) != null){
							Integer existingCount = this.wordCount.get(arr[3]);
							Integer newCount = existingCount + Integer.parseInt(arr[0]);
							this.wordCount.put(arr[3],newCount);
						}
						else{
							this.wordCount.put(word, Integer.parseInt(arr[0]));
						}
					}
					//setting the number of times a word has a tag sequence of "O"
					if(arr[2] != null && arr[2].trim().length() > 0 && arr[2].toLowerCase().equals("o")){
						if(oMap.get(arr[3]) != null){
							String oWord = arr[3];
							Integer existingOCount = oMap.get(oWord);
							Integer newOCount = existingOCount + Integer.parseInt(arr[0]);
							oMap.put(oWord,newOCount);
						}
						if(oMap.get(arr[3]) == null){
							oMap.put(arr[3], Integer.parseInt(arr[0]));
						}
					}
					//setting the number of times a word has a tag sequence of "I-GENE"
					if(arr[2] != null && arr[2].trim().length() > 0 && arr[2].toLowerCase().equals("i-gene")){
						if(iGeneMap.get(arr[3]) != null){
							String iGeneWord = arr[3];
							Integer existingIGeneCount = iGeneMap.get(iGeneWord);
							Integer newIGeneCount = existingIGeneCount + Integer.parseInt(arr[0]);
							iGeneMap.put(iGeneWord,newIGeneCount);
						}
						if(iGeneMap.get(arr[3]) == null){
							iGeneMap.put(arr[3], Integer.parseInt(arr[0]));
						}
					}					
				}
				// Setting the unigram count
				if(arr[1].toLowerCase().equals("1-gram")){
					if(arr[2].toLowerCase().equals("i-gene")){
						this.iGeneSequenceCount = Integer.parseInt(arr[0]);
					}
					else{
						this.oSequenceCount = Integer.parseInt(arr[0]);
					}
				}
				// SEtting the bigram count
				if(arr[1].toLowerCase().equals("2-gram")){
					String key = arr[2].toLowerCase()+","+arr[3].toLowerCase();
					Integer value = Integer.parseInt(arr[0]);
					this.bigramCount.put(key, value);
				}
				//Setting the trigram count
				if(arr[1].toLowerCase().equals("3-gram")){
					String key = arr[2].toLowerCase()+","+arr[3].toLowerCase()+","+arr[4].toLowerCase();
					Integer value = Integer.parseInt(arr[0]);
					this.triGramCount.put(key, value);
				}	
			}
			bufferReader.close();
			fileReader.close();
		}
		catch(FileNotFoundException ex){
			System.out.println(ex);	
		}
		catch(IOException ex){
			System.out.println("error reading file" + fileName);
		}
	}
	
	/* Reads the training or testing file and regenerate gene.counts file after replacing rare words and 
	 words that have numbers, allcaps and lastcap  */
	private void replaceRareWords(){
		String fileName = "gene.train";
		String line = null;
		try {
			FileReader fileReader = new FileReader(fileName);
			BufferedReader bufferReader = new BufferedReader(fileReader);
			List<String> lines = new ArrayList<String>();
			Path file = Paths.get("temp.train");
			while((line = bufferReader.readLine()) != null){
				if(line.trim().length() != 0){
					String arr[] = line.split(" ");
					String newLine = null;
					boolean allCaps = true;
					boolean hasNumber = false;
					boolean lastCap = false;
					
					for (int i = 0; i < arr[0].length(); i++){
					    char c = arr[0].charAt(i);        
					    int asciiVal= (int) c;
					    if(asciiVal > 47 && asciiVal < 58){
					    	hasNumber = true;
					    } else if(asciiVal < 65 || asciiVal > 90){
					    	 allCaps = false;
						} else if(i == arr[0].length() -1){
					    	char last = arr[0].charAt(i);
					    	int asciiLastVal= (int) last;
					    	if(asciiLastVal >= 65 && asciiLastVal <= 90){
					    		lastCap = true;
							} 	
					    }
					}
					
					if((this.wordCount.get(arr[0]) == null || this.wordCount.get(arr[0]) < 5) && hasNumber){
						arr[0] ="_HASNUMBER_";
						newLine = arr[0] +" "+arr[1];
						lines.add(newLine);
					}else if((this.wordCount.get(arr[0]) == null || this.wordCount.get(arr[0]) < 5) && allCaps){
						arr[0] ="_ALLCAPS_";
						newLine = arr[0] +" "+arr[1];
						lines.add(newLine);
					} else if((this.wordCount.get(arr[0]) == null || this.wordCount.get(arr[0]) < 5) && lastCap){
						arr[0] ="_LASTCAP_";
						newLine = arr[0] +" "+arr[1];
						lines.add(newLine);
					} else if((this.wordCount.get(arr[0]) == null || this.wordCount.get(arr[0]) < 5)){
						arr[0]="_RARE_";
						newLine = arr[0] +" "+arr[1];
						lines.add(newLine);
					} else {
						newLine = arr[0] +" "+arr[1];
						lines.add(newLine);
					}
				}
				else{
					String newLine="";
					lines.add(newLine);
				}
			}
			Path oldTrainFile = Paths.get("gene.train");
			Files.delete(oldTrainFile);
			Files.write(file, lines, Charset.forName("UTF-8"));
			new File("temp.train").renameTo(new File("gene.train"));
			Path oldGeneFile = Paths.get("gene.counts");
			Files.delete(oldGeneFile);
			String[] cmd = {
			        "/bin/bash",
			        "-c",
			        "python count_freqs.py gene.train > gene.counts",
		    };
			try {
				Runtime.getRuntime().exec(cmd).waitFor();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			this.mapWords();
			bufferReader.close();
			fileReader.close();
		}
		catch(FileNotFoundException ex){
			System.out.println(ex);
		}
		catch(IOException ex){
			System.out.println(ex);
		}
	}
	
	public Integer getOSequenceCount() {
		return oSequenceCount;
	}

	public void setOSequenceCount(Integer oSequenceCount) {
		this.oSequenceCount = oSequenceCount;
	}

	public Integer getiGeneSequenceCount() {
		return iGeneSequenceCount;
	}

	public void setiGeneSequenceCount(Integer iGeneSequenceCount) {
		this.iGeneSequenceCount = iGeneSequenceCount;
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
