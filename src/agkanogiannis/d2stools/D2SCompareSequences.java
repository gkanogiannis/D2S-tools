package agkanogiannis.d2stools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.tree.DefaultTreeModel;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.io.IOUtils;

import com.google.common.collect.ImmutableList;

import agkanogiannis.Coder;
import agkanogiannis.Utils;
import agkanogiannis.hcluster.Clade;
import agkanogiannis.hcluster.HierarchicalCluster;
import agkanogiannis.hcluster.PhylipWriter;
import agkanogiannis.io.FastaManager;

public class D2SCompareSequences {
	
	private static String version;
	
	private FastaManager frm;
	
	public static void main(String[] args) {
		D2SCompareSequences cmp = new D2SCompareSequences();
		version = new Date(Utils.classBuildTimeMillis(cmp.getClass())).toString();
		System.out.println("version D2SCompareSequences ="+version);
		
		int numOfThreads = 1;
    	int k = 4;
    	List<String> inputFileNames = null;
    	String outputPrefix = null;
    	boolean normalize = false;
    	
    	ArrayList<String> seqNames = new ArrayList<String>();
    	ArrayList<ReadD2> seqVectors = new ArrayList<ReadD2>();
    	
	    CommandLineParser parser = new BasicParser();
	    Options options = cmp.createOptions();
	    HelpFormatter formatter = new HelpFormatter();
	    CommandLine cmd = null;
	    try {
	        cmd = parser.parse( options, args );
	        
	        if(args.length<1){
	    		throw new Exception();
	    	}

	        if(cmd.hasOption("t")){
	        	numOfThreads = Integer.parseInt(cmd.getOptionValue("t"));
	        }
	        if(cmd.hasOption("k")){
	        	k = Integer.parseInt(cmd.getOptionValue("k"));
	        }
	        if(cmd.hasOption("i")){
	        	inputFileNames = Arrays.asList(cmd.getOptionValues("i"));
	        }
	        if(cmd.hasOption("o")){
	        	outputPrefix = cmd.getOptionValue("o");
	        }
	        if(cmd.hasOption("n")){
	        	if(cmd.getOptionValue("n").equalsIgnoreCase("yes"))
	        		normalize = true;
	        	else if(cmd.getOptionValue("n").equalsIgnoreCase("no"))
	        		normalize = false;
	        }
	        
	        
	        for(String inputFileName : inputFileNames){
		    	cmp.createVectorsForSequences(numOfThreads, k, inputFileName, seqNames, seqVectors, normalize);
		    }
		    
		    double[][] distances = cmp.calculateComparisonsSequences(seqNames, seqVectors);
		    
		    cmp.saveMatrix(seqNames, distances, outputPrefix);
		    
		    cmp.hclustering(seqNames.toArray(new String[seqNames.size()]), distances, outputPrefix, 10);
	    }
	    catch(Exception exp ) {
	    	System.out.println(exp.getMessage());
	    	formatter.printHelp( "D2SCompareSequences", options );
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
		Option input = OptionBuilder.withArgName("input").withLongOpt("input").hasArgs(Integer.MAX_VALUE).withDescription("Input Fasta/q files paths.").isRequired(true).create("i");
		options.addOption(input);
		Option output = OptionBuilder.withArgName("output").withLongOpt("output").hasArg().withDescription("Output files prefix.").isRequired(true).create("o");
		options.addOption(output);
		Option normalize = OptionBuilder.withArgName("normalize").withLongOpt("normalize").hasArg().withDescription("Normalize vectors, yes|no").isRequired(false).create("n");
		options.addOption(normalize);
		
		return options;
	}
	
	private void createVectorsForSequences(int numOfThreads, int k, String inputFileName, ArrayList<String> seqNames, ArrayList<ReadD2> seqVectors, boolean normalize){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(1+1);
		    //CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
		    
		    //System.out.println(Utils.time()+" START of Counting for Sequence="+seqName);
		    
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads+1);
			
			if(frm!=null){
				frm.clear();
				frm = null;
			}
			frm = new FastaManager(false, ImmutableList.of(inputFileName), startSignal, doneSignal);
			pool.execute(frm);
			
			
			startSignal.await();
			boolean done = false;
			AtomicInteger seqCount = new AtomicInteger(0);
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
				seqCount.incrementAndGet();
				//String seqName = new String(read.getHeader()).split("\\s++")[0].substring(1);
				String seqName = new String(read.getHeader()).substring(1).
							         replaceAll("[^A-Za-z0-9]", " ").
							         replaceAll(" +", " ");
				System.out.println(seqName);
				if(seqName.isEmpty() || seqName==null ) seqName = String.valueOf(read.getReadId());
				seqNames.add(seqName);
				
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
				//calculate per sequence probs
				read.calculateProbs(k);
				if(normalize) {
					read.normalizeProbs(read.getNorm());
				}
				ReadD2 seqVector = read;
				//ReadD2Centroid seqVector = new ReadD2Centroid(read);
				seqVectors.add(seqVector);
				//read.clear();
								
				//System.out.println("No. reads="+seqVector.getNumOfElements());
				System.out.println("No. Counts="+seqVector.getTotalCounts());
				System.out.println("No. ATCG="+seqVector.getTotalATCG());
				System.out.println("As="+seqVector.getAs());
				System.out.println("Ts="+seqVector.getTs());
				System.out.println("Cs="+seqVector.getCs());
				System.out.println("Gs="+seqVector.getGs());
				System.out.println("norm="+seqVector.getNorm());
				System.out.println(Utils.time()+" END of Vector Creation for Sequence="+seqName);
				System.out.println();
			}
			doneSignal.countDown();
			doneSignal.await();
			pool.shutdown();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private double[][] calculateComparisonsSequences(List<String> seqNames, List<ReadD2> seqVectors){
		try{			
			int numOfSeqToCompare = seqNames.size();
			
			double[][] distances = new double[numOfSeqToCompare][numOfSeqToCompare];
			
			for(int i=0; i<numOfSeqToCompare; i++){
				ReadD2 X = seqVectors.get(i);
				for(int j=i; j<numOfSeqToCompare; j++){
					ReadD2 Y = seqVectors.get(j);
					double d2_measure = DissimilarityMeasuresD2.d2_S_Dissimilarity(X, Y);
					if(Math.abs(d2_measure) < Double.valueOf("1E-15")) d2_measure=0.0;
					distances[i][j] = d2_measure;
					distances[j][i] = d2_measure;
					//System.out.println(Utils.time()+"\td2S for couple ["+seqNames.get(i)+" :: "+seqNames.get(j)+"]="+d2_measure);
				}
			}
			
			return distances;
		}
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private void saveMatrix(List<String> seqNames, double[][] distances, String outputPrefix){
		try{
			File f = new File(outputPrefix+".d2S.dist").getCanonicalFile();
			f.getParentFile().mkdirs();
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));
			bw.write(seqNames.size()+"\n");
			for(int i=0; i<seqNames.size(); i++){
				bw.write(seqNames.get(i).replaceAll(" +", "_"));
				for(int j=0; j<seqNames.size(); j++){
					bw.write("\t"+distances[i][j]);
				}
				bw.write("\n");
			}
			
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void hclustering(String[] seqNames, double[][] distances, String outputPrefix, int minClusterSize){
		try {
			System.out.println(Utils.time()+" Distances="+distances.length+"x"+distances[0].length);
			
			String method = HierarchicalCluster.AVERAGE;
			HierarchicalCluster hc = new HierarchicalCluster(distances, seqNames);
			hc.setLinkageMethod(method);
			Clade root = hc.cluster();
			
			PhylipWriter writer = new PhylipWriter();
			StringWriter sw = new StringWriter();
			writer.setOutput(sw);
			try {
				writer.write(new DefaultTreeModel(root));
			} 
			catch(IOException ex) {
				ex.printStackTrace();
			}
			
			String treeString = sw.toString().replace("\n", "");
			System.out.println("hierarchical method="+method);
			FileOutputStream fos = new FileOutputStream(new File(outputPrefix+".tree"));
			IOUtils.write(treeString+"\n", fos);
			fos.close();
			
			String[] labelsReordered = Utils.reorderLabels(seqNames, treeString);
			double[][] distancesReordered = Utils.reorderDistances(distances, seqNames, labelsReordered);
		
			System.out.println("preJRI");
			JRITools jritools = JRITools.getInstance(null);
			System.out.println("precut");
			TreeMap<Integer, TreeSet<String>> clusters = jritools.dynamicTreeCut(treeString, distancesReordered, labelsReordered, minClusterSize);
			System.out.println("afterCut");
			jritools.shutdown();
			
			File f = new File(outputPrefix+".d2S-hclusters.txt").getCanonicalFile();
			f.getParentFile().mkdirs();
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));
			System.out.println(Utils.time()+" Clusters="+clusters.size()+"\n");
			bw.write("Clusters="+clusters.size()+"\n");
			for(Entry<Integer, TreeSet<String>> entry : clusters.entrySet()){
				int clusterId = entry.getKey();
				TreeSet<String> cluster = entry.getValue();
				System.out.println("Cluster "+clusterId+"="+cluster.size());
				bw.write("Cluster "+clusterId+"="+cluster.size()+"\n");
				for(String name : cluster){
					System.out.println("\t"+name);
					bw.write("\t"+name+"\n");
				}
			}
			bw.flush();
			bw.close();
		}
		catch(Exception e) {
		    e.printStackTrace();
		}
	}
	
}
