package agkanogiannis.d2stools;

import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.hash.TLongIntHashMap;

import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import agkanogiannis.Utils;

public class Read {

	private final int readId;
	
	private byte[] header;
	private byte[] seq;
	private byte[] qual;

	private int abundanceCluster = -1;
	
	private TLongIntHashMap kmerCounts;
	
	private int[] ranks = null;
	
	public Read(int readId, byte[] header, byte[] seq) {
		this.readId = readId;
		this.header = header;
		this.seq = seq;
		kmerCounts = new TLongIntHashMap();
	}
	
	public Read(int readId, byte[] header, byte[] seq, byte[] qual) {
		this.readId = readId;
		this.header = header;
		this.seq = seq;
		this.qual = qual;
		kmerCounts = new TLongIntHashMap();
	}
	
	public void initRanks(Map<Long, Integer> spaceRanks){
		if(ranks == null){
			int spaceSize = spaceRanks.size();
			ranks = new int[spaceSize];
		
			//put present kmers in a hashmap
			TreeMap<Long, Integer> unsorted = new TreeMap<Long, Integer>();
			for ( TLongIntIterator it = iteratorCounts(); it.hasNext(); ) {
				it.advance();
				long kmerCode = it.key();
				int kmerCount = it.value();
				unsorted.put(kmerCode, kmerCount);
			}
			//System.out.println("unsorted1 "+unsorted.size());
			//System.out.println(unsorted);
			
			//insert zeros for non present kmers
			for(long kmerCode : spaceRanks.keySet()){
				if(!unsorted.containsKey(kmerCode)){
					unsorted.put(kmerCode, 0);
				}
			}
			//System.out.println("unsorted2 "+unsorted.size());
			//System.out.println(unsorted);
			
			//sort the map
			Map<Long, Integer> sorted =  Utils.sortByValue(unsorted);	
			//System.out.println("sorted "+sorted.size());
			//System.out.println(sorted);
			
			//write ranks
			int current = 0;
			for(Entry<Long, Integer> entry : sorted.entrySet()){
				ranks[spaceRanks.get(entry.getKey())] = current++;
			}
			//System.out.println("ranks "+ranks.length);
			//System.out.println(Arrays.toString(ranks));
			//System.out.println("space "+spaceRanks.size());
			//System.out.println(spaceRanks.toString());
			
			unsorted.clear();
			sorted.clear();
			unsorted = null;
			sorted = null;
			//System.exit(0);
		}
	}
	
	public int[] getRanks(){
		return ranks;
	}
	
	/*
	public void remove(ArrayList<Long> removed){
		for(long remove : removed){
			kmers.remove(remove);
		}
	}
	*/
	
	public void insertKmerCount(long kmerCode, int count){
		if(kmerCode<0L){
			return;
		}
		kmerCounts.adjustOrPutValue(kmerCode, count, count);
	}
	
	public int getCountForKmerCode(long kmerCode){
		return kmerCounts.get(kmerCode);
	}
	
	public void clearHead(){
		header = null;
	}
	
	public void clearSeq(){
		seq = null;
		qual = null;
	}
	
	public void clearHeadSeq(){
		header = null;
		seq = null;
		qual = null;
	}

	public void clearVec(){
		if(kmerCounts!=null){
			kmerCounts.clear();
			kmerCounts = null;
		}
	}
	
	public void clearFull(){
		header = null;
		seq = null;
		qual = null;
		clearVec();
	}
	
	public int getReadId() {
		return readId;
	}

	public byte[] getHeader(){
		return header;
	}
	
	public byte[] getSeq(){
		return seq;
	}
	
	public byte[] getQual(){
		return qual;
	}
	
	public int getAbundanceCluster() {
		return abundanceCluster;
	}

	public void setAbundanceCluster(int abundanceCluster) {
		this.abundanceCluster = abundanceCluster;
	}

	public char charAt(int index) {
		return (char) (seq[index] & 0xff);
	}

	public int getLength() {
		return seq==null?0:seq.length;
	}
	
	public void addWith(Read other){
		for ( TLongIntIterator it = other.iteratorCounts(); it.hasNext(); ) {
			it.advance();
			insertKmerCount(it.key(), it.value());
		}
	}
	
	public TLongIntIterator iteratorCounts() {
		return new MyIterator();
	}
	
	private class MyIterator implements TLongIntIterator {

        private final TLongIntIterator iterator;
        private int current;
        private boolean hasNext;

        MyIterator() {
        	iterator = kmerCounts.iterator();
        	current = kmerCounts.size();
        }

        public boolean hasNext() {
        	hasNext = (current > 0);
            return hasNext;
        }

        public void advance() {
			if(!hasNext){
    			return;
    		}
        	iterator.advance();
            current--;
		}
        
        public long key() {
        	return iterator.key();
		}
        
        public int value() {
        	return iterator.value();
		}
        
        
        //Unused
        public void remove() {
      	}
      	public int setValue(int arg0) {
      		return 0;
      	}	
      	
	}
	
}
