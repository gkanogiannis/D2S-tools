/*
 *
 * D2S-tools agkanogiannis.d2stools.D2SCompareReads
 *
 * Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package agkanogiannis.d2stools;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.io.IOUtils;

import agkanogiannis.Utils;
import agkanogiannis.io.FastaManager;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

public class D2SCompareReads {
	
	private static String version;
	
	//for long kmer filtering
	private static DictionaryFilter dictionaryFilter;
	private HashMap<Integer, ReadProcessorD2> readProcessors;
	private Map<Integer, ReadD2> readsD2;
	
	public static void main(String[] args) {
		D2SCompareReads cmp = new D2SCompareReads();
		version = new Date(Utils.classBuildTimeMillis(cmp.getClass())).toString();
		System.out.println("version D2SCompareReads ="+version);
		
		int numOfThreads = 1;
    	int k = 4;
    	//kmer length for filtering
    	int kFilter = 20;
    	List<String> inputFileNames = null;
    	String outputFileName = null;
    	
	    CommandLineParser parser = new BasicParser();
	    Options options = cmp.createOptions();
	    HelpFormatter formatter = new HelpFormatter();
	    CommandLine cmd = null;
	    try {
	        cmd = parser.parse( options, args );

	        if(cmd.hasOption("t")){
	        	numOfThreads = Integer.parseInt(cmd.getOptionValue("t"));
	        }
	        if(cmd.hasOption("k")){
	        	k = Integer.parseInt(cmd.getOptionValue("k"));
	        }
	        if(cmd.hasOption("kF")){
	        	kFilter = Integer.parseInt(cmd.getOptionValue("kF"));
	        	if(kFilter<=1){
	        		kFilter = 1;
	        	}
	        }
	        if(cmd.hasOption("i")){
	        	inputFileNames = Arrays.asList(cmd.getOptionValues("i"));
	        }
	        if(cmd.hasOption("o")){
	        	outputFileName = cmd.getOptionValue("o");
	        } 
	        
	        if(inputFileNames.size() > 1){
		    	cmp.processReads(numOfThreads, k, kFilter, inputFileNames);
			    //cmp.createReadVectors(numOfThreads);
			    cmp.calculateAndSaveComparisonsReadsReads(numOfThreads, outputFileName);
		    }
		    else if(inputFileNames.size() == 1) {
		    	
		    }
		    else{
	    		throw new Exception("Please give at least one input file.");
	    	}
	    }
	    catch(Exception exp ) {
	    	System.out.println(exp.getMessage());
	    	formatter.printHelp( "D2SCompareReads", options );
	    	System.exit(0);
	    }
	    
	    System.out.println(Utils.RAMInfo(Runtime.getRuntime()));
	}
	
	@SuppressWarnings("static-access")
	private Options createOptions(){
		Options options = new Options();
		Option t = OptionBuilder.withArgName("numOfThreads").withLongOpt("numOfThreads").hasArg().withDescription("Number of threads to use.").isRequired(false).create("t");
		options.addOption(t);
		Option k = OptionBuilder.withArgName("kMerSize").withLongOpt("kMerSize").hasArg().withDescription("k-mer length.").isRequired(false).create("k");
		options.addOption(k);
		Option kFilter = OptionBuilder.withArgName("kMerSizeFilter").withLongOpt("kMerSizeFilter").hasArg().withDescription("k-mer length for filtering.").isRequired(false).create("kF");
		options.addOption(kFilter);
		Option input = OptionBuilder.withArgName("input").withLongOpt("input").hasArgs(Integer.MAX_VALUE).withDescription("Input Reads Fasta/q files paths.").isRequired(true).create("i");
		options.addOption(input);
		Option output_D = OptionBuilder.withArgName("output").withLongOpt("output").hasArg().withDescription("Output files prefix.").isRequired(true).create("o");
		options.addOption(output_D);
		return options;
	}
	
	private void processReads(int numOfThreads, int k, int kFilter, List<String> inputFileNames){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
		    dictionaryFilter = new DictionaryFilter(1024*1024, usingThreads);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
		    
		    System.out.println(Utils.time()+" START of Counting Reads D2");
		    
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads+1);
			
			readProcessors = new  HashMap<Integer, ReadProcessorD2>();
			FastaManager frm = new FastaManager(false, inputFileNames, startSignal, doneSignal);
			pool.execute(frm);
			
			ReadProcessorD2.resetCounters();
			//for d2S filtering
			//int kk = 20;
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				ReadProcessorD2 rp = new ReadProcessorD2(null, dictionaryFilter, ReadProcessorD2.MODE.KMER_COUNTING_READ_D2, k, kFilter, frm, startSignal, doneSignal);
				readProcessors.put(rp.getId(), rp);
				pool.execute(rp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" END of Counting Reads D2");
			System.out.println(Utils.time()+" Loaded reads: "+ReadProcessorD2.getReadCount().get());
			
			System.out.println(Utils.time()+" START populate read2read relation");
			dictionaryFilter.populateRead2ReadRelation(usingThreads);
			System.out.println(Utils.time()+" END populate read2read relation");
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
    /*
	private void createReadVectors(int numOfThreads){
		try{
			System.out.println(Utils.time()+" START of Vector Creation for Reads");
			
			readsD2 = new HashMap<Integer, ReadD2>();
			//long As = 0;
			//long Ts = 0;
			//long Cs = 0;
			//long Gs = 0;
			for(ReadProcessorD2 rp : readProcessors.values()){
				//As += rp.As;
				//Ts += rp.Ts;
				//Cs += rp.Cs;
				//Gs += rp.Gs;
				for(ReadD2 readD2 : rp.getReadsD2()){
					readsD2.put(readD2.getReadId(), readD2);
				}
			}
				
			
			//Parallel for
			System.out.println("PARALLEL loop");
    		//final long _As = As;
    		//final long _Ts = Ts;
    		//final long _Cs = Cs;
    		//final long _Gs = Gs;

			for(ReadProcessorD2 rp : readProcessors.values()){
			    int usingThreads = (Runtime.getRuntime().availableProcessors() < numOfThreads ? Runtime.getRuntime().availableProcessors() : numOfThreads);
			    ParallelFor.blockingFor(usingThreads, rp.getReadsD2(), 
			    		 // The operation to perform with each item
			    		 new ParallelFor.Operation<ReadD2>() {
			    		    public void perform(ReadD2 readD2) {
			    		    	readD2.calculateProbs(0, 0, 0, 0);
			    		    	//readD2.calculateProbs(_As, _Ts, _Cs, _Gs);
			    		    	//readD2.normalize();
								//System.out.println(readD2.getReadId()+"\t"+readD2.getNorm());
			    		    };
			    		});
			}
			
		    
			System.out.println("NORMAL loop");
			for(ReadProcessorD2 rp : readProcessors.values()){
				for(ReadD2 readD2 : rp.getReadsD2()){
					readD2.calculateProbs(As, Ts, Cs, Gs);
					readsD2.add(readD2);
					//System.out.println("No. Length="+readD2.getLength());
					//System.out.println("No. ATCG="+readD2.getTotalATCG());
					//System.out.println("No. ATCG2="+readD2.getTotalATCG(k));
				}
			}
			
			
			System.out.println(Utils.time()+" END of Vector Creation for Reads");
			System.out.println(Utils.time()+" Proccessed reads: "+readsD2.size());
			System.out.println();	
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
    */
	
	private void calculateAndSaveComparisonsReadsReads(int numOfThreads, String outputFileName){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads);
			
			//Starting threads
			ArrayList<ReadReadD2CalculatorThread> threads = new ArrayList<ReadReadD2CalculatorThread>();
			int numOfReads = readsD2.size();
			double chunk = (double)numOfReads / (double)usingThreads;
			System.out.println("chunk="+chunk);
			int is = 1;
			for(int i=0; i<usingThreads; i++){
				int ie = (int)Math.floor((double)is+chunk);
				if(i==usingThreads-1 || ie>numOfReads){
					ie = numOfReads;
				}
				System.out.println("Thread "+(i+1)+" "+(is)+"->"+(ie));
				ReadReadD2CalculatorThread t = new ReadReadD2CalculatorThread(i+1, is, ie, doneSignal);
				threads.add(t);
				pool.execute(t);
				is = ie+1;
			}
			doneSignal.await();
			pool.shutdown();
			
			System.out.println("Writing file.");
			DecimalFormat df = new DecimalFormat("#.#############");
			File f = new File(outputFileName).getCanonicalFile();
	 		f.getParentFile().mkdirs();
	 		BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(f));
	 		for(ReadReadD2CalculatorThread t : threads){
				for(int i=0; i<t.arrayRead1.size(); i++){
					IOUtils.write(t.arrayRead1.get(i)+"\t"+t.arrayRead2.get(i)+"\t"+df.format(1.0-t.arrayMeasure.get(i))+"\n", bos);
				}
			}
	 		bos.flush(); bos.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private class ReadReadD2CalculatorThread implements Runnable {
		private int id;
		private int istart, iend;
		private CountDownLatch doneSignal;
		
		TIntArrayList arrayRead1 = new TIntArrayList();
		TIntArrayList arrayRead2 = new TIntArrayList();
		TDoubleArrayList arrayMeasure = new TDoubleArrayList();
		
		private ReadReadD2CalculatorThread(int id, int istart, int iend, CountDownLatch doneSignal){
			this.id = id;
			this.istart = istart;
			this.iend = iend;
			this.doneSignal = doneSignal;
		}
		
		@Override
		public void run() {
			int sum = iend-istart;
			//System.out.println(sum);
			ReadD2 read1, read2;
			double d2_measure;
			
			int currPercentage = 0;
			int percentage;
			int currSum = 0;
			
			for(int i=istart;i<=iend;i++){
				read1 = readsD2.get(i);
				for(int read2ID : dictionaryFilter.getR2RListForId(read1.getReadId())){
					read2 = readsD2.get(read2ID);
					//if(dictionaryFilter.areRelated(read1, read2)){
						d2_measure = DissimilarityMeasuresD2.d2_S_Dissimilarity(read1, read2);
						arrayRead1.add(read1.getReadId());
						arrayRead2.add(read2.getReadId());
						arrayMeasure.add(d2_measure);
					//}
				}
				currSum++;
				percentage = (int)(100.0*((double)currSum/(double)sum));
				if(percentage % 25 == 0 && percentage > currPercentage){
					currPercentage = percentage;
					System.out.println(Utils.time()+" Thread "+id+" "+percentage+"%");
				}
			}
			doneSignal.countDown();
		}
	}
	
}
