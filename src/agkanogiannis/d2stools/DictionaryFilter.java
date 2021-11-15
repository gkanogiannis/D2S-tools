package agkanogiannis.d2stools;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentMap;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.MapMaker;

import agkanogiannis.Utils;

public class DictionaryFilter {

	private  ConcurrentMap<Long, TIntList> kmer2readMap;
	private  HashMultimap<Integer, Integer> read2readMap;
	
	public DictionaryFilter(int initialCapacity, int concurrencyLevel){
		kmer2readMap = new MapMaker().concurrencyLevel(concurrencyLevel).initialCapacity(initialCapacity).makeMap();
	}
	
	/*
	public HashMultimap<Integer, Integer> getR2R(){
		return read2readMap;
	}
	*/
	
	public List<Integer> getR2RListForId(int readId){
		return Utils.asSortedList(read2readMap.get(readId));
	}
	
	public Set<Integer> getR2RSetForId(int readId){
		return read2readMap.get(readId);
	}
	
	public synchronized void insert(long kmerId, final int readId){
		if(kmerId<0L){
			return;
		}
		
		//Collection<Integer> old = kmer2readMap.putIfAbsent(kmerId, new ConcurrentSkipListSet<Integer>(ImmutableList.of(readId)));
		//Collection<Integer> old = kmer2readMap.putIfAbsent(kmerId, new ConcurrentLinkedDeque<Integer>(ImmutableList.of(readId)));
		TIntList old = kmer2readMap.putIfAbsent(kmerId, new TIntArrayList(){{add(readId);}});
		if(old != null) {
			if(!old.contains(readId)){
				old.add(readId);
			}
		}
		
		return;
	}
	
	public void populateRead2ReadRelation(int usingThreads){
		if(read2readMap!=null){
			read2readMap.clear();
			read2readMap = null;
		}
		read2readMap = HashMultimap.create();
			
		/*
		//PARALLEL
		System.out.println("PARALLEL");
		ParallelForTrove.blockingFor(usingThreads, kmer2readMap.values(), 
	    		 // The operation to perform with each item
	    		 new ParallelForTrove.Operation<TIntList>() {
	    		    public void perform(TIntList list) {
	    		    	list.sort();
	    				for(int i=0; i<list.size(); i++){
	    					for(int j=i; j<list.size(); j++){
	    						read2readMap.get(list.get(i)).add(list.get(j));
	    					}
	    				}
	    		    };
	    		});
		*/
		
		
		//SERIAL
		for(TIntList list : kmer2readMap.values()){
			//System.out.println("sortinglist "+(++count)+" size="+list.size());
			//list.sort();
			for(int i=0; i<list.size(); i++){
				for(int j=i+1; j<list.size(); j++){
					read2readMap.get(list.get(i)).add(list.get(j));
				}
			}
		}
	}
	
	public boolean areRelated(Read read1, Read read2){
		return read2readMap.get(read1.getReadId()).contains(read2.getReadId()) || read2readMap.get(read2.getReadId()).contains(read1.getReadId());
	}
	
	public boolean areRelated(int read1Id, int read2Id){
		return read2readMap.get(read1Id).contains(read2Id) || read2readMap.get(read2Id).contains(read1Id);
	}
	
	public String toString(int kFilter) {
		StringBuilder sb = new StringBuilder();
		int numOflists = kmer2readMap.values().size();
		sb.append("numOflists:\t\t"+numOflists+"\n");
		
		return sb.toString();
	}
	
}
