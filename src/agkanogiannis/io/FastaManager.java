package agkanogiannis.io;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedDeque;
//import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.CountDownLatch;

import agkanogiannis.Utils;
import agkanogiannis.d2stools.Read;

public class FastaManager implements Runnable{
	private ArrayList<Read> reads = null;
	
	private ConcurrentLinkedDeque<Integer> readsIds = null;
	//private ConcurrentSkipListSet<Integer> readsIds = null;
	
	private int numOfReads = 0;
	
	public boolean isFastq = false;
	private boolean keepQualities = false;
	private List<String> inputFileNames = null;
	private boolean done = false;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	//private BTreeMap<Integer, Read> reads = null;
	
	public FastaManager(boolean keepQualities, List<String> inputFileNames, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.keepQualities = keepQualities;
		this.inputFileNames = inputFileNames;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		
		this.reads = new ArrayList<Read>();
		
		this.readsIds = new ConcurrentLinkedDeque<Integer>();
		//this.readsIds = new ConcurrentSkipListSet<Integer>();
		
		//this.reads = Utils.getDB().getTreeMap("reads");// getHashMap("reads");
	}
	
	public void clear(){
		if(reads!=null){
			reads.clear();
		}
		reads = null;
		
		if(readsIds!=null){
			readsIds.clear();
		}
		readsIds = null;
	}
	
	public List<Read> getReads() {
		return Collections.unmodifiableList(reads);
	}
	
	public boolean hasMore() {
		if(!done){
			return true;
		}
		else{
			if(this.readsIds.isEmpty()){
				return false;
			}
			else{
				return true;
			}
		}
	}
	
	private void putRead(Read read) {
		/*
		System.out.println("id="+read.getReadId());
		System.out.println("length="+read.getLength());
		if(read.getHeader()!=null)
		System.out.println("header="+new String(read.getHeader()));
		if(read.getSeq()!=null)
		System.out.println(" seq="+new String(read.getSeq()));
		if(read.getQual()!=null)
		System.out.println("qual="+new String(read.getQual()));
		//System.exit(0);
		*/
		
		reads.add(read);
		readsIds.add(read.getReadId());
	}
	
	public Read getConcurrentRead() {
		Integer readId = readsIds.pollFirst();
		if(readId == null){
			return null;
		}
		else{
			return reads.get(readId-1);
		}
	}
	
	public void run() {
		try{
			done = false;
			
		    System.out.println(Utils.time()+" FastaManager: START READ");
		    startSignal.countDown();
		    
		    List<File> inputFiles = new ArrayList<File>();
		    for(String inputFileName : inputFileNames){
		    	inputFileName = inputFileName.trim();
		    	if( !isFastq && (inputFileName.endsWith(".fastq") || inputFileName.endsWith(".fq"))){
		    		isFastq = true;
		    	}
		    	inputFiles.add(new File(inputFileName).getCanonicalFile());
		    }
		    
		    //Check if files exist
		    for(File f : inputFiles){
		    	if(!f.exists() || !f.canRead()){
		    		System.out.println("\tERROR : File "+f+"\n\tdoes not exist ot cannot be read. Exiting.");
		    		System.exit(1);
		    	}
		    }
		    
		    FastaDecoder decoder = new FastaDecoder();
		    FastaIterator<byte[]> iterator = FastaIterator.create(decoder, inputFiles);
		    byte[] lastHeader = null;
		    ArrayList<byte[]> lastSeq = null;
		    byte[] lastQualSeq = null;
		    for (List<byte[]> chunk : iterator) {
		    	System.out.println(Utils.time()+" FastaManager: Chunk with "+chunk.size()+" lines.");
		    	int i = 0;
		    	
		    	//From previous chunk
		    	if(lastHeader!=null){
		    		
		    		//FASTA
		    		if(!isFastq){
		    			while( i < chunk.size()){
		    				byte[] line = chunk.get(i);
		    				//Read sequence
							if(line[0] != '>'){
								lastSeq.add(line);
							}
							else{
								break;
							}
							
							i++;
		    			}
		    			//System.out.println(numOfReads+" From previous chunk A");
		    			//putRead(new Read(numOfReads++, lastHeader, Utils.concat(lastSeq)));
		    			putRead(new Read(++numOfReads, lastHeader, Utils.concat(lastSeq)));
		    			lastHeader = null;
		    		}
		    		
		    		//FASTQ
		    		else{
		    			byte[] seq = null;
		    			byte[] qual = null;
						if(lastQualSeq!=null){
							seq = lastQualSeq;
							
							if(chunk.get(i)[0]=='+'){
								i++;
							}
							//Read qualities
							qual = chunk.get(i);
						}
						else{
							//Read sequence
							seq = chunk.get(i);
								
							i++;
							i++;
							//Read qualities
							qual = chunk.get(i);
						}
						//System.out.println(numOfReads+" From previous chunk Q");
						if(keepQualities){
			    			//putRead(new Read(numOfReads++, lastHeader, seq, qual));
							putRead(new Read(++numOfReads, lastHeader, seq, qual));
						}
			    		else{
			    			//putRead(new Read(numOfReads++, lastHeader, seq));
			    			putRead(new Read(++numOfReads, lastHeader, seq));
			    		}
						lastHeader = null;
						lastQualSeq = null;
						
						i++;
		    		}	
		    		
		    	}
		    	
		    	
		    	while( i < chunk.size() ){
		    		byte[] line = chunk.get(i);
		    		//System.out.println(i+" "+chunk.size()+" "+new String(line));
		    		
		    		//FASTA
		    		if(!isFastq && line[0] == '>'){
		    			ArrayList<byte[]> seq = new ArrayList<byte[]>();
		    			byte[] header = line;
						while( ++i < chunk.size()){
							line = chunk.get(i);
							//Read sequence
							if(line[0] != '>'){
								seq.add(line);
							}
							else{
								break;
							}
						}
						
						if(i == chunk.size()){
							lastHeader = header;
							lastSeq = seq;
						}
						else{
							//System.out.println(numOfReads+" normalA");
							//putRead(new Read(numOfReads++, header, Utils.concat(seq)));
							putRead(new Read(++numOfReads, header, Utils.concat(seq)));
						}
					}
		    			
		    		//FASTQ
		    		else if(isFastq && line[0] == '@'){
		    			byte[] seq = null;
		    			byte[] qual = null;
		    			byte[] header = line;
						if( ++i < chunk.size()){
							//Read sequence
							seq = chunk.get(i);
							
							if( ++i < chunk.size()){
								if( ++i < chunk.size()){
									//Read qualities
									qual = chunk.get(i);
								}
							}
						}
						
						if(i == chunk.size()){
							lastHeader = header;
							lastQualSeq = seq;
						}
						else{
							//System.out.println(numOfReads+" normalQ");
				    		if(keepQualities){
				    			//putRead(new Read(numOfReads++, header, seq, qual));
				    			putRead(new Read(++numOfReads, header, seq, qual));
				    		}
				    		else{
				    			//putRead(new Read(numOfReads++, header, seq));
				    			putRead(new Read(++numOfReads, header, seq));
				    		}
						}
						
						i++;
					}
		    		
		    		else{
		    			i++;
		    		}
		    	
		    	}
		    	chunk.clear();
		    	chunk = null;
		    }
		    
		    
		    //From last chunk
	    	if(lastHeader!=null){
	 		
	    		//FASTA
	    		if(!isFastq){	
	    			//System.out.println(numOfReads+" From last chunk A");
	    			//putRead(new Read(numOfReads++, lastHeader, Utils.concat(lastSeq)));
	    			putRead(new Read(++numOfReads, lastHeader, Utils.concat(lastSeq)));
	    			lastHeader = null;
	    		}
	    		
	    		//FASTQ
	    		else{
					if(!keepQualities){
						//System.out.println(numOfReads+" From last chunk Q");
		    			//putRead(new Read(numOfReads++, lastHeader, lastQualSeq));
						putRead(new Read(++numOfReads, lastHeader, lastQualSeq));
					}
					lastHeader = null;
					lastQualSeq = null;
	    		}	
	    	}
		    
			System.out.println(Utils.time()+" FastaManager: END READ");
			System.out.println(Utils.time()+" FastaManager: "+(isFastq?"FASTQ":"FASTA"));
			done = true;
			//startSignal.countDown();
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

}
