/*
 *
 * D2S-tools agkanogiannis.d2stools.JRITools
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

import java.awt.FileDialog;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.TreeSet;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

public class JRITools {
	private static JRITools instance = null;
	private static Rengine re = null;

	@SuppressWarnings("unused")
	private JRITools() {
	}

	public static JRITools getInstance(String[] args) {
		if (instance == null) {
			instance = new JRITools(args);
		}
		return instance;
	}

	public JRITools(String[] args) {
		try {
			if(args==null){
				args = new String[]{"--vanilla"};
			}
			// just making sure we have the right version of everything
			if (!Rengine.versionCheck()) {
				System.err.println("** Version mismatch - Java files don't match library version.");
				System.exit(1);
			}

			System.out.println("Creating Rengine (with arguments ) = "+Arrays.toString(args));
			re = new Rengine(args, false, new TextConsole());
			System.out.println("Rengine created, waiting for R");
			// the engine creates R is a new thread, so we should wait until
			// it's ready
			if (!re.waitForR()) {
				System.out.println("Cannot load R");
				return;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void shutdown() {
		try {
			re.end();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public TreeMap<Integer, TreeSet<String>> dynamicTreeCut(String treeString, double[][] distances, String[] labels, int minClusterSize) {
		try {
			REXP x;
			re.eval("library(\"dynamicTreeCut\")");
			re.eval("library(\"ape\")");
			//re.eval("library(\"dendextend\")");

			re.assign("treeString", treeString);
			//System.out.println("Got back="+re.eval("treeString").asString());
			
			assignAsRMatrix(re, distances, "distancesMatrix");
			re.assign("rawLabels", labels);
			re.assign("colLabels", labels);
			re.eval("rownames(distancesMatrix) <- rawLabels");
			re.eval("colnames(distancesMatrix) <- colLabels");
			re.eval("distancesAsDist <- as.dist(distancesMatrix)");
			//System.out.println("rawLabels="+re.eval("rawLabels"));
			//System.out.println("colLabels="+re.eval("colLabels"));
			//System.out.println("\ndistancesMatrix=\n"+re.eval("print(distancesMatrix)"));
			//System.out.println("\ndistancesAsDist=\n"+re.eval("print(distancesAsDist)"));
			//System.out.println(treeString);
			
			re.eval("tree <- read.tree(text=treeString);");
			System.out.println("\nIs ultrametric=" + re.eval("is.ultrametric(tree);").asBool().toString());

			//re.eval("tree.ultra <- chronos(tree);");
			//System.out.println("\nIs ultrametric=" + re.eval("is.ultrametric(tree.ultra);").asBool().toString());
			
			re.eval("tree_hclust <- as.hclust(tree);");
			
			/*
			String command = "treecut <- cutreeDynamic(dendro=tree_hclust, method=\"tree\", cutHeight=NULL, minClusterSize=3, deepSplit=TRUE);";
			re.eval(command);
			System.out.println("\ntreecut");
			System.out.println(x=re.eval("treecut"));
			double[] result = x.asDoubleArray();
			*/
			
			String cutHeight = "0.005";
			
			//String command = "treecut_hybrid <- cutreeDynamic(dendro=tree_hclust, method=\"hybrid\", distM=distancesMatrix, cutHeight=NULL, minClusterSize=1, deepSplit=TRUE, maxCoreScatter=NULL, minGap=NULL, maxAbsCoreScatter=NULL, minAbsGap=NULL);";
			String command = "treecut_hybrid <- cutreeDynamic(dendro = tree_hclust, cutHeight = "+cutHeight+", minClusterSize = " + minClusterSize + ", method = \"hybrid\", distM = distancesMatrix, deepSplit = 1, maxCoreScatter = NULL, minGap = NULL, " + 
							 "pamStage = TRUE, pamRespectsDendro = TRUE, maxDistToLabel = NULL, respectSmallClusters = TRUE, " + 
							 "verbose = 2, indent = 0);";
			x=re.eval(command);
			System.out.println("\ntreecut_hybrid");
			System.out.println(x=re.eval("treecut_hybrid"));
			double[] result = x.asDoubleArray();
			
			
			/*
			x = re.eval("tree_hclust$order");
			int[] ordersFromTreeCut = x.asIntArray();
			x = re.eval("tree_hclust$labels");
			String[] labelsFromTreeCut = x.asStringArray();
			for (int i = 0; i < ordersFromTreeCut.length; i++) {
				System.out.println("order="+ordersFromTreeCut[i]+" "+(labels[i].equalsIgnoreCase(labelsFromTreeCut[i])?"OK":"ERROR")+" "+" label="+labels[i]+" labelFromTreeCut="+labelsFromTreeCut[i]);
			}
			*/
			
			
			//find final clusters
			TreeMap<Integer, TreeSet<String>> clusters = new TreeMap<Integer, TreeSet<String>>();
			for(int i=0; i<result.length; i++){
				int clusterId = (int)result[i];
				TreeSet<String> cluster = clusters.get(clusterId);
				if(cluster==null){
					cluster = new TreeSet<String>();
					clusters.put(clusterId, cluster);
				}
				String label = labels[i];
				cluster.add(label);
			}
			return clusters;
			
			//return null;
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	private REXP assignAsRMatrix(Rengine rEngine, double[][] sourceArray, String nameToAssignOn) {
        if (sourceArray.length == 0) {
            return null;
        } 
        
        /*
        rEngine.assign(nameToAssignOn, sourceArray[0]);
        REXP resultMatrix = rEngine.eval(nameToAssignOn + " <- matrix(" + nameToAssignOn + " ,nr=1)");
        for (int i = 1; i < sourceArray.length; i++) {
            rEngine.assign("temp", sourceArray[i]);
            resultMatrix = rEngine.eval(nameToAssignOn + " <- rbind(" + nameToAssignOn + ",matrix(temp,nr=1))");
            System.out.println(i);
        }
        */
        
        //flatten
        double[] temp1Darray = new double[sourceArray.length * sourceArray[0].length];
        for(int i=0; i<sourceArray.length; i++){
        	for(int j=0; j<sourceArray[i].length; j++){
        		temp1Darray[i+j] = sourceArray[i][j];
        	}
        }
        
        rEngine.assign("temp", temp1Darray);
        REXP resultMatrix = rEngine.eval(nameToAssignOn + " <- matrix(temp, " + sourceArray.length + ", " + sourceArray[0].length + ", byrow = T) ");
       
        return resultMatrix;
    }

	private class TextConsole implements RMainLoopCallbacks {
		public void rWriteConsole(Rengine re, String text, int oType) {
			System.out.print(text);
		}

		public void rBusy(Rengine re, int which) {
			System.out.println("rBusy(" + which + ")");
		}

		public String rReadConsole(Rengine re, String prompt, int addToHistory) {
			System.out.print(prompt);
			try {
				BufferedReader br = new BufferedReader(new InputStreamReader(
						System.in));
				String s = br.readLine();
				return (s == null || s.length() == 0) ? s : s + "\n";
			} catch (Exception e) {
				System.out.println("jriReadConsole exception: " + e.getMessage());
			}
			return null;
		}

		public void rShowMessage(Rengine re, String message) {
			System.out.println("rShowMessage \"" + message + "\"");
		}

		@SuppressWarnings("deprecation")
		public String rChooseFile(Rengine re, int newFile) {
			FileDialog fd = new FileDialog(new Frame(),
					(newFile == 0) ? "Select a file" : "Select a new file",
					(newFile == 0) ? FileDialog.LOAD : FileDialog.SAVE);
			fd.show();
			String res = null;
			if (fd.getDirectory() != null)
				res = fd.getDirectory();
			if (fd.getFile() != null)
				res = (res == null) ? fd.getFile() : (res + fd.getFile());
			return res;
		}

		public void rFlushConsole(Rengine re) {
		}

		public void rLoadHistory(Rengine re, String filename) {
		}

		public void rSaveHistory(Rengine re, String filename) {
		}
	}
	
}
