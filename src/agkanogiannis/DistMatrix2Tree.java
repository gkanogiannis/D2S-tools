/*
 *
 * D2S-tools agkanogiannis.DistMatrix2Tree
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
package agkanogiannis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringWriter;

import javax.swing.tree.DefaultTreeModel;

import agkanogiannis.hcluster.Clade;
import agkanogiannis.hcluster.HierarchicalCluster;
import agkanogiannis.hcluster.PhylipWriter;

public class DistMatrix2Tree {

	public static void main(String[] args) {
		DistMatrix2Tree cmd = new DistMatrix2Tree();
		try {
			BufferedReader br = new BufferedReader(new FileReader(args[0]));
			String line = br.readLine();
			System.err.println(line);
			int numOfSamples = Integer.parseInt(line.split("[\\s,\\t]+")[0]);
			double[][] distances = new double[numOfSamples][numOfSamples];
			String[] sampleNames = new String[numOfSamples];
			int i = 0;
			while((line=br.readLine())!=null) {
				String[] data = line.split("[\\s,\\t]+");
				sampleNames[i] = data[0].trim();
				for(int j=0; j<numOfSamples; j++) {
					try {
						distances[i][j] = Double.parseDouble(data[j+1]);
					} 
					catch (NumberFormatException e) {
						distances[i][j] = 1.0;
					}
				}
				i++;
			}
			br.close();
			
			cmd.hclustering(sampleNames, distances);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void hclustering(String[] sampleNames, double[][] distances){
		try {
			System.err.println(Utils.time()+" Distances="+distances.length+"x"+distances[0].length);
			
			String method = HierarchicalCluster.COMPLETE;
			HierarchicalCluster hc = new HierarchicalCluster(distances, sampleNames);
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
			System.err.println("hierarchical method="+method);
			System.out.println(treeString);
		}
		catch(Exception e) {
		    e.printStackTrace();
		}
	}

}
