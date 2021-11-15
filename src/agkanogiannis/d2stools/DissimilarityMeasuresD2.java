package agkanogiannis.d2stools;

import gnu.trove.iterator.TLongIntIterator;

import java.util.HashSet;

public class DissimilarityMeasuresD2 {
	
	public static double d2_S_Dissimilarity(ReadD2Interface X, ReadD2Interface Y){
		double D2_S = 0.0;
		double tempX = 0.0;
		double tempY = 0.0;
		double cXi = 0.0;
		double cYi = 0.0;
		double cXi_bar = 0.0;
		double cYi_bar = 0.0;
		double temp3 = 0.0;
		
		HashSet<Long> set = new HashSet<Long>();
		
		for ( TLongIntIterator itX = X.iteratorCounts(); itX.hasNext(); ) {
			itX.advance();
			long kmerCodeX = itX.key();
			set.add(kmerCodeX);
		}
		for ( TLongIntIterator itY = Y.iteratorCounts(); itY.hasNext(); ) {
			itY.advance();
			long kmerCodeY = itY.key();
			set.add(kmerCodeY);
		}
		//set.addAll(X.getKmerCodes());
		//set.addAll(Y.getKmerCodes());
		
		//Iterator<Long> iter = X.iterator();
		//long kmerCode;
		for(long kmerCode : set){
		//while(iter.hasNext()){
			//kmerCode = iter.next();
			/*
			cXi = X.getCountForKmerCode(kmerCode);
			cYi = Y.getCountForKmerCode(kmerCode);
			if(cYi > 0){
				pXi = X.getProbForKmerCode(kmerCode);
				pYi = Y.getProbForKmerCode(kmerCode);
				cXi_bar = (double)cXi - (double)(X.getLength()-k+1)*pXi;
				cYi_bar = (double)cYi - (double)(Y.getLength()-k+1)*pYi;
				temp1 = Math.pow(cXi_bar, 2.0);
				temp2 = Math.pow(cYi_bar, 2.0);
				temp3 = Math.sqrt(temp1+temp2);
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				D2_S += cXi_bar*cYi_bar/temp3;
				tempX += cXi_bar*cXi_bar/temp3;
				tempY += cYi_bar*cYi_bar/temp3;
			}
			*/
			cYi = Y.getDoubleCountForKmerCode(kmerCode);
			cXi = X.getDoubleCountForKmerCode(kmerCode);
			//System.out.println("kmerCode="+kmerCode);
			//System.out.println("\tcountProbX="+X.getDoubleCountForKmerCode(kmerCode)+"  :  "+X.getDoubleProbForKmerCode(kmerCode));
			//System.out.println("\tcountProbY="+Y.getDoubleCountForKmerCode(kmerCode)+"  :  "+Y.getDoubleProbForKmerCode(kmerCode));
			//if(cYi > 0.0){
				cXi_bar = cXi - (double)X.getTotalCounts()*X.getDoubleProbForKmerCode(kmerCode);
				cYi_bar = cYi - (double)Y.getTotalCounts()*Y.getDoubleProbForKmerCode(kmerCode);
				//cXi_bar = X.getCountForKmerCodeDouble(kmerCode) - X.getProbForKmerCode(kmerCode);
				//cYi_bar = cYi - Y.getProbForKmerCode(kmerCode);
				//System.out.println("\tcXi_bar="+ cXi_bar);
				//System.out.println("\tcYi_bar="+ cYi_bar);
				temp3 = Math.sqrt(Math.pow(cXi_bar, 2.0) + Math.pow(cYi_bar, 2.0));
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				D2_S += (cXi_bar*cYi_bar)/temp3;
				tempX += (cXi_bar*cXi_bar)/temp3;
				tempY += (cYi_bar*cYi_bar)/temp3;
			//}
		}
		
		//System.out.println(tempY);
		
		tempX = Math.sqrt(tempX);
		tempY = Math.sqrt(tempY);
		double temp = D2_S/(tempX*tempY);
		return 0.5*(1.0 - temp);
	}
	
	public static double d2_Star_Dissimilarity(ReadD2Interface X, ReadD2Interface Y){
		double D2_Star = 0.0;
		double tempX = 0.0;
		double tempY = 0.0;
		double cYi = 0;
		double cXi_bar = 0.0;
		double cYi_bar = 0.0;
		double temp3 = 0.0;
		
		HashSet<Long> set = new HashSet<Long>();
		
		for ( TLongIntIterator itX = X.iteratorCounts(); itX.hasNext(); ) {
			itX.advance();
			long kmerCodeX = itX.key();
			set.add(kmerCodeX);
		}
		for ( TLongIntIterator itY = Y.iteratorCounts(); itY.hasNext(); ) {
			itY.advance();
			long kmerCodeY = itY.key();
			set.add(kmerCodeY);
		}
		//set.addAll(X.getKmerCodes());
		//set.addAll(Y.getKmerCodes());
		
		for(long kmerCode : set){
			cYi = Y.getDoubleCountForKmerCode(kmerCode);
			//if(cYi > 0.0){
				cXi_bar = X.getDoubleCountForKmerCode(kmerCode) - (double)X.getTotalCounts()*X.getDoubleProbForKmerCode(kmerCode);
				cYi_bar = cYi - (double)Y.getTotalCounts()*Y.getDoubleProbForKmerCode(kmerCode);
				temp3 = Math.sqrt((double)X.getTotalCounts()*X.getDoubleProbForKmerCode(kmerCode)) * Math.sqrt((double)Y.getTotalCounts()*Y.getDoubleProbForKmerCode(kmerCode));
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				D2_Star += (cXi_bar*cYi_bar)/temp3;
				tempX += (cXi_bar*cXi_bar)/temp3;
				tempY += (cYi_bar*cYi_bar)/temp3;
			//}
		}
		
		tempX = Math.sqrt(tempX);
		tempY = Math.sqrt(tempY);
		double temp = D2_Star/(tempX*tempY);
		return 0.5*(1.0 - temp);
	}
	
	/*
	public static double d2_Star_Dissimilarity(ReadD2 X, ReadD2 Y, int k){
		double D2_Star = 0.0;
		double d2_Star = 0.0;
		double tempX = 0.0;
		double tempY = 0.0;
		int cXi = 0;
		int cYi = 0;
		double pXi = 0.0;
		double pYi = 0.0;
		double cXi_bar = 0.0;
		double cYi_bar = 0.0;
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		
		Iterator<Long> iter = X.iterator();
		long kmerCode;
		while(iter.hasNext()){
			kmerCode = iter.next();
			cXi = X.getCountForKmerCode(kmerCode);
			cYi = Y.getCountForKmerCode(kmerCode);
			if(cYi > 0){
				pXi = X.getProbForKmerCode(kmerCode);
				pYi = Y.getProbForKmerCode(kmerCode);
				temp1 = (double)(X.getLength()-k+1)*pXi;
				temp2 = (double)(Y.getLength()-k+1)*pYi;
				cXi_bar = (double)cXi - temp1;
				cYi_bar = (double)cYi - temp2;
				temp3 = Math.sqrt(temp1*temp2);
				if(temp1 == 0.0){
					temp1 = 1.0;
					temp3 = 1.0;
				}
				if(temp2 == 0.0){
					temp2 = 1.0;
					temp3 = 1.0;
				}
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				
				D2_Star += cXi_bar*cYi_bar/temp3;
				tempX += cXi_bar*cXi_bar/temp1;
				tempY += cYi_bar*cYi_bar/temp2;
			}
		}
		
		tempX = Math.sqrt(tempX);
		tempY = Math.sqrt(tempY);
		double temp = D2_Star/(tempX*tempY);
		d2_Star = 0.5*(1-temp);
		return d2_Star;
	}
	*/

}
