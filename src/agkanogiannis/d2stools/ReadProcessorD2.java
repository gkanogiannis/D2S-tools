package agkanogiannis.d2stools;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import agkanogiannis.Coder;
import agkanogiannis.Utils;
import agkanogiannis.io.FastaManager;

public class ReadProcessorD2 implements Runnable {

	private static AtomicInteger readCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private Dictionary dictionary = null;
	//For filtering with long kFilter
	private DictionaryFilter dictionaryFilter = null;
	private MODE mode = null;
	private int k;
	//For filtering with long kFilter
	private int kFilter = 0;
	private FastaManager frm = null;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	private ReadD2Centroid rpSampleVector = null;
	
	public static enum MODE {
		KMER_COUNTING_READ_D2,
		KMER_COUNTING_SAMPLE_D2
	}
	
	public ReadD2Centroid getrpSampleVector(){
		return rpSampleVector;
	}
	
	public static void resetCounters(){
		readCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public ReadProcessorD2(Dictionary dictionary, DictionaryFilter dictionaryFilter, MODE mode, int k, int kFilter, FastaManager frm, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.dictionary = dictionary;
		this.dictionaryFilter = dictionaryFilter;
		this.mode = mode;
		this.k = k;
		this.kFilter = kFilter;
		this.frm = frm;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
	}
	
	public static AtomicInteger getReadCount() {
		return readCount;
	}

	public MODE getMode() {
		return mode;
	}

	public void setStartSignal(CountDownLatch startSignal) {
		this.startSignal = startSignal;
	}

	public void setDoneSignal(CountDownLatch doneSignal) {
		this.doneSignal = doneSignal;
	}

	public void setMode(MODE mode) {
		this.mode = mode;
	}

	public int getId() {
		return id;
	}

	public void run() {
		try{
			if(rpSampleVector==null){
				rpSampleVector = new ReadD2Centroid(new Read(id, null, null));
			}
			
			startSignal.await();
			boolean done = false;
			//System.out.println(Utils.time()+" ReadProcessorD2: "+id+"\t"+mode.toString()+" START");
			switch(mode){
				case KMER_COUNTING_READ_D2 :
					while(!done){
						Read r = frm.getConcurrentRead();
						ReadD2 read = null;
						if(r!=null){
							read = new ReadD2(r);
						}
						if(read==null){
							if(!frm.hasMore()){
								done = true;
								break;
							}
							continue;
						}
						processRead_KMER_COUNTING_READ_D2(read);
					}
					break;
				case KMER_COUNTING_SAMPLE_D2 :
					while(!done){
						Read r = frm.getConcurrentRead();
						ReadD2 read = null;
						if(r!=null){
							read = new ReadD2(r);
						}
						if(read==null){
							if(!frm.hasMore()){
								done = true;
								break;
							}
							continue;
						}
						processRead_KMER_COUNTING_SAMPLE_D2(read);
					}
					break;
				default:
					System.err.println(Utils.time()+" Mode+\"" + mode.toString() + "\" not known.");
					break;
			}
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
	private void processRead_KMER_COUNTING_READ_D2(ReadD2 read){
		try{
			readCount.incrementAndGet(); 
			if((read.getReadId()) % 1000000 == 0){
				System.out.println(Utils.time()+" ReadProcessorD2: "+id+"\t"+mode.toString()+" "+(read.getReadId()));
			}
			
			long kmerCode;
			long kmerCodeFilter;	
			
			for(int i=0; i+k<=read.getLength(); i++){	
				int oldAs = read.as;
				int oldTs = read.ts;
				int oldCs = read.cs;
				int oldGs = read.gs;
				kmerCode = Coder.encodeToLong(read, i, i+k, false, true);
				if(kmerCode>=0L){
					if(dictionary!=null){
						dictionary.insert(kmerCode);
					}
					read.insertKmerCount(kmerCode, 1);
					read.insertKmerProb(kmerCode, (short)(read.as-oldAs), (short)(read.ts-oldTs), (short)(read.cs-oldCs), (short)(read.gs-oldGs));
				}
				
				//reverse
				oldAs = read.as;
				oldTs = read.ts;
				oldCs = read.cs;
				oldGs = read.gs;
				kmerCode = Coder.encodeToLong(read, i, i+k, true, true);
				if(kmerCode>=0L){
					if(dictionary!=null){
						dictionary.insert(kmerCode);
					}
					read.insertKmerCount(kmerCode, 1);
					read.insertKmerProb(kmerCode, (short)(read.as-oldAs), (short)(read.ts-oldTs), (short)(read.cs-oldCs), (short)(read.gs-oldGs));
				}
			}
		
			//calculate per read probs
			read.calculateProbs(k);
			rpSampleVector.addWith((ReadD2Interface)read, true);
			read.clear();
			

			/////////////////////
			//long kmer filtering
			if(kFilter>0 && dictionaryFilter!=null){
				for(int i=0; i+kFilter<=read.getLength(); i++){	
					kmerCodeFilter = Coder.encodeToLong(read, i, i+kFilter, false, false);
					if(kmerCodeFilter>=0L){
						dictionaryFilter.insert(kmerCodeFilter, read.getReadId());
					}
					
					//reverse
					kmerCodeFilter = Coder.encodeToLong(read, i, i+kFilter, true, false);
					if(kmerCodeFilter>=0L){
						dictionaryFilter.insert(kmerCodeFilter, read.getReadId());
					}
				}
			}
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	private void processRead_KMER_COUNTING_SAMPLE_D2(ReadD2 read){
		try{
			readCount.incrementAndGet();
			if((read.getReadId()) % 1000000 == 0){
				System.out.println(Utils.time()+" ReadProcessorCMG: "+id+"\t"+mode.toString()+" "+(read.getReadId()));
			}
			
			long kmerCode;
			
			for(int i=0; i+k<=read.getLength(); i++){
				int oldAs = read.as;
				int oldTs = read.ts;
				int oldCs = read.cs;
				int oldGs = read.gs;
				kmerCode = Coder.encodeToLong(read, i, i+k, false, true);
				if(kmerCode>=0L){
					read.insertKmerCount(kmerCode, 1);
					read.insertKmerProb(kmerCode, (short)(read.as-oldAs), (short)(read.ts-oldTs), (short)(read.cs-oldCs), (short)(read.gs-oldGs));
				}
				
				//reverse
				oldAs = read.as;
				oldTs = read.ts;
				oldCs = read.cs;
				oldGs = read.gs;
				kmerCode = Coder.encodeToLong(read, i, i+k, true, true);
				if(kmerCode>=0L){
					read.insertKmerCount(kmerCode, 1);
					read.insertKmerProb(kmerCode, (short)(read.as-oldAs), (short)(read.ts-oldTs), (short)(read.cs-oldCs), (short)(read.gs-oldGs));
				}
			}
			
			rpSampleVector.addWith((ReadD2Interface)read, false);
			read.clear();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
}
