package agkanogiannis.d2stools;

import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.iterator.TLongIntIterator;

public interface ReadD2Interface {
	
	public double  getDoubleCountForKmerCode(long kmerCode);
	public double  getDoubleProbForKmerCode(long kmerCode);
	
	public TLongIntIterator    iteratorCounts();
	public TLongDoubleIterator iteratorProbs();
	
	public long getTotalATCG();
	public long getTotalCounts();
	
	public long getNumOfElements();
	
	public long getAs();
	public long getTs();
	public long getCs();
	public long getGs();
}