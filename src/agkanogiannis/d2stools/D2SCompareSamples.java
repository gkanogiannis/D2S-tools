package agkanogiannis.d2stools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.tree.DefaultTreeModel;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;

import agkanogiannis.Utils;
import agkanogiannis.hcluster.Clade;
import agkanogiannis.hcluster.HierarchicalCluster;
import agkanogiannis.hcluster.PhylipWriter;
import agkanogiannis.io.FastaManager;

public class D2SCompareSamples {
	
	private static String version;
	
	private FastaManager frm;
	private HashMap<Integer, ReadProcessorD2> readProcessors;
	
	public static void main(String[] args) {
		D2SCompareSamples cmp = new D2SCompareSamples();
		version = new Date(Utils.classBuildTimeMillis(cmp.getClass())).toString();
		System.out.println("version D2SCompareSamples ="+version);
		
		int numOfThreads = 1;
    	int k = 4;
    	List<String> inputFileNames = null;
    	String outputPrefix = null;
    	
    	ArrayList<String> samplesNames = new ArrayList<String>();
    	ArrayList<ArrayList<String>> samplesFiles = new ArrayList<ArrayList<String>>();
    	ArrayList<ReadD2Centroid> samplesVectors = new ArrayList<ReadD2Centroid>();
    	
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
	        
	        if(inputFileNames.size()==1){
	        	// Sample files info in txt file with format 
	        	// Name1:path_set1.1.fq.gz;path_set1.2.fq.gz;path_set1.3.fq;...
		    	if(!inputFileNames.get(0).endsWith(".fa") && !inputFileNames.get(0).endsWith(".fa.gz") &&
		    	   !inputFileNames.get(0).endsWith(".fasta") && !inputFileNames.get(0).endsWith(".fasta.gz") &&
		    	   !inputFileNames.get(0).endsWith(".fq") && !inputFileNames.get(0).endsWith(".fq.gz") &&
		    	   !inputFileNames.get(0).endsWith(".fastq") && !inputFileNames.get(0).endsWith(".fastq.gz") )
		    	{
		    		BufferedReader br = new BufferedReader(new FileReader(inputFileNames.get(0)));
		    		String line = null;
		    		//Name1:path_set1.1.fq.gz;path_set1.2.fq.gz;path_set1.3.fq;...
		    		while((line=br.readLine())!=null){
		    			ArrayList<String> list = new ArrayList<String>();
		    			samplesFiles.add(list);
		    			samplesNames.add(line.split(":")[0].trim());
		    			String[] files = line.split(":")[1].trim().split(";");
		    			for(String file : files){
		    				list.add(file);
		    			}
		    		}
		    		br.close();
		    	}
		    	else{
		    		throw new Exception("Please give more than one input files.");
		    	}
		    }
	        else{
	        	for(int i=1; i<=inputFileNames.size(); i++){
	        		samplesNames.add(FilenameUtils.removeExtension(inputFileNames.get(i-1)));
	        		ArrayList<String> list = new ArrayList<String>();
	        		list.add(inputFileNames.get(i-1));
	        		samplesFiles.add(list);
	        	}
	        }
	        
	        for(int i=0; i<samplesNames.size(); i++){
		    	samplesVectors.add(cmp.createVectorForSample(numOfThreads, k, samplesFiles.get(i), samplesNames.get(i), (i+1)));
		    }
		    
		    double[][] distances = cmp.calculateComparisonsSamples(samplesNames, samplesVectors);
		    
		    cmp.saveMatrix(samplesNames, distances, outputPrefix);
		    
		    cmp.hclustering(samplesNames.toArray(new String[samplesNames.size()]), distances, outputPrefix, 10);
	    }
	    catch(Exception exp ) {
	    	System.out.println(exp.getMessage());
	    	formatter.printHelp( "D2SCompareSamples", options );
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
		
		return options;
	}
	
	private ReadD2Centroid createVectorForSample(int numOfThreads, int k, ArrayList<String> sampleFiles, String sampleName, int sampleId){
		try{
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.out.println("cpus="+cpus);
		    System.out.println("using="+usingThreads);
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
		    
		    System.out.println(Utils.time()+" START of Counting for Sample="+sampleName);
		    
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads+1);
			
			readProcessors = new  HashMap<Integer, ReadProcessorD2>();
			if(frm!=null){
				frm.clear();
				frm = null;
			}
			frm = new FastaManager(false, sampleFiles, startSignal, doneSignal);
			pool.execute(frm);
			
			ReadProcessorD2.resetCounters();
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				ReadProcessorD2 rp = new ReadProcessorD2(null, null, ReadProcessorD2.MODE.KMER_COUNTING_SAMPLE_D2, k, 0, frm, startSignal, doneSignal);
				readProcessors.put(rp.getId(), rp);
				pool.execute(rp);
			}
			
			doneSignal.await();
			pool.shutdown();
			
			System.out.println(Utils.time()+" END of Counting for Sample="+sampleName);
			System.out.println(Utils.time()+" Loaded reads: "+ReadProcessorD2.getReadCount().get());
			
			
			System.out.println(Utils.time()+" START of Vector Creation for Sample="+sampleName);
			ReadD2Centroid sampleVector = new ReadD2Centroid(new Read(-sampleId, null, null));
	
			for(ReadProcessorD2 rp : readProcessors.values()){
				sampleVector.addWith((ReadD2Interface)rp.getrpSampleVector(), false);
			}
		
			sampleVector.calculateProbs(k);
			
			System.out.println("No. reads="+sampleVector.getNumOfElements());
			System.out.println("No. Counts="+sampleVector.getTotalCounts());
			System.out.println("No. ATCG="+sampleVector.getTotalATCG());
			System.out.println("As="+sampleVector.getAs());
			System.out.println("Ts="+sampleVector.getTs());
			System.out.println("Cs="+sampleVector.getCs());
			System.out.println("Gs="+sampleVector.getGs());
			//System.out.println("norm="+sampleVector.getNorm());
			System.out.println(Utils.time()+" END of Vector Creation for Sample="+sampleName);
			System.out.println();
			
			return sampleVector;
		}
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private double[][] calculateComparisonsSamples(List<String> samplesNames, List<ReadD2Centroid> samplesVectors){
		try{			
			int numOfSamplesToCompare = samplesNames.size();
			
			double[][] distances = new double[numOfSamplesToCompare][numOfSamplesToCompare];
			
			for(int i=0; i<numOfSamplesToCompare; i++){
				ReadD2Centroid X = samplesVectors.get(i);
				for(int j=i; j<numOfSamplesToCompare; j++){
					ReadD2Centroid Y = samplesVectors.get(j);
					double d2_measure = DissimilarityMeasuresD2.d2_S_Dissimilarity(X, Y);
					if(Math.abs(d2_measure) < Double.valueOf("1E-15")) d2_measure=0.0;
					distances[i][j] = d2_measure;
					distances[j][i] = d2_measure;
					System.out.println(Utils.time()+"\td2S for couple ["+samplesNames.get(i)+":"+samplesNames.get(j)+"]="+d2_measure);
				}
			}
			
			return distances;
		}
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private void saveMatrix(List<String> samplesNames, double[][] distances, String outputPrefix){
		try{
			File f = new File(outputPrefix+".d2S-matrix.csv").getCanonicalFile();
			f.getParentFile().mkdirs();
			BufferedWriter bw = new BufferedWriter(new FileWriter(f));
			bw.write("Sample");
			for(int i=0; i<samplesNames.size(); i++){
				bw.write(","+samplesNames.get(i));
			}
			bw.write("\n");
			for(int i=0; i<samplesNames.size(); i++){
				for(int j=0; j<samplesNames.size(); j++){
					if(j==0){
						bw.write(samplesNames.get(i));
					}
					bw.write(","+(distances[i][j]));
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
	
	private void hclustering(String[] samplesNames, double[][] distances, String outputPrefix, int minClusterSize){
		try {
			System.out.println(Utils.time()+" Distances="+distances.length+"x"+distances[0].length);
			
			String method = HierarchicalCluster.AVERAGE;
			HierarchicalCluster hc = new HierarchicalCluster(distances, samplesNames);
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
			
			String[] labelsReordered = Utils.reorderLabels(samplesNames, treeString);
			double[][] distancesReordered = Utils.reorderDistances(distances, samplesNames, labelsReordered);
		
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
