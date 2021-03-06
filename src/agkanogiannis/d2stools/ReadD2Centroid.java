/*
 *
 * D2S-tools agkanogiannis.d2stools.ReadD2Centroid
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

import java.nio.ByteBuffer;

import gnu.trove.iterator.TLongDoubleIterator;

public class ReadD2Centroid extends ReadD2 implements ReadD2Interface{

	private long numOfElements = 0L;
	
	private long As = 0L;
	private long Ts = 0L;
	private long Cs = 0L;
	private long Gs = 0L;
	
	public ReadD2Centroid(Read read) {
		super(read);
	}
	
	public long getAs() {return As;}
	public long getTs() {return Ts;}
	public long getCs() {return Cs;}
	public long getGs() {return Gs;}

	public long getNumOfElements() {
		return numOfElements;
	}
	
	public long getTotalATCG() {
		return As+Ts+Cs+Gs;
	}
	
	public void addWith(ReadD2Interface other, boolean adjust){
		if(other==null){
			return;
		}
		if(other instanceof ReadD2Centroid){
			//System.out.println("Adding with numOfElements="+other.getNumOfElements());
		}
		super.addWith((ReadD2)other);
		numOfElements += other.getNumOfElements();
		As += other.getAs();
		Ts += other.getTs();
		Cs += other.getCs();
		Gs += other.getGs();
		for ( TLongDoubleIterator it = other.iteratorProbs(); it.hasNext(); ) {
			it.advance();
			if(adjust){
				adjustKmerProb(it.key(), it.value());
			}
			else{
				insertKmerProb(it.key(), it.value());
			}
		}
	}
	
	public void calculateProbs(int k){
		double temp = 0.0;
		long Total = getTotalATCG();
		//System.out.println("id:"+getReadId()+"\tTotal="+Total);
		ByteBuffer bb = ByteBuffer.wrap(new byte[8]);
		for ( TLongDoubleIterator it = iteratorProbs(); it.hasNext(); ) {
			it.advance();
			long kmerCode = it.key();
			
			((ByteBuffer) bb.position(0)).putDouble(getDoubleProbForKmerCode(kmerCode)).position(0);
			
			short a = bb.getShort();
			short t = bb.getShort();
			short c = bb.getShort();
			short g = bb.getShort();
				
			double prob  =  Math.pow((double)As/(double)Total, (double)a);
			prob *=  		Math.pow((double)Ts/(double)Total, (double)t);
			prob *=  		Math.pow((double)Cs/(double)Total, (double)c);
			prob *=  		Math.pow((double)Gs/(double)Total, (double)g);
				
			insertKmerProbForce(kmerCode, prob);
			temp += prob;
		}
		System.out.println("sumprob="+temp);
		//System.out.println("distinct kmers="+kmerProbs.size());
	}

}
