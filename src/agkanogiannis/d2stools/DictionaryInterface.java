package agkanogiannis.d2stools;

import agkanogiannis.hcluster.ClusterPoisson;
import agkanogiannis.hcluster.ClusterVectorTrove;
import gnu.trove.map.hash.TIntLongHashMap;

public interface DictionaryInterface {

	public int getMaxCount();
	public TIntLongHashMap getCountsHisto();
	public void insert(long kmerCode);
	public int getGlobalCountFor(long kmerCode);
	public void clear();
	public void removeAll(int excludeMin);
	public ClusterVectorTrove[] createABClusterVectors(ClusterPoisson[] clusterPoissons, int excludeMin, int excludeMax);
}
