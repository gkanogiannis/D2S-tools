/*
 *
 * D2S-tools agkanogiannis.d2stools.DictionaryInterface
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
