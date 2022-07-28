package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.Callable;

import population.Population_Bridging;
import relationship.ContactMap;

public abstract class Abstract_Runnable_ClusterModel implements Runnable {

	protected String runnableId = null;
	protected File baseDir = null;

	public String getRunnableId() {
		return runnableId;
	}

	public void setRunnableId(String runnableId) {
		this.runnableId = runnableId;
	}

	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	public abstract Object[] getRunnable_fields();

	public static final Callable<ContactMap> generateContactMapCallable(File cMap_file) {
		Callable<ContactMap> callable = new Callable<ContactMap>() {
			@Override
			public ContactMap call() throws Exception {
				StringWriter cMap_str = new StringWriter();
				PrintWriter pWri = new PrintWriter(cMap_str);
				BufferedReader reader = new BufferedReader(new FileReader(cMap_file));
				String line;
				while ((line = reader.readLine()) != null) {
					pWri.println(line);
				}

				pWri.close();
				reader.close();
				return ContactMap.ContactMapFromFullString(cMap_str.toString());
			}

		};
		return callable;
	}

	public static final Callable<ContactMap> generateContactMapCallable(String cMap_str) {
		Callable<ContactMap> callable = new Callable<ContactMap>() {
			@Override
			public ContactMap call() throws Exception {
				return ContactMap.ContactMapFromFullString(cMap_str);
			}
		};
		return callable;
	}

	public static final Callable<ArrayList<Integer[]>> generateContactEdgeArray(ContactMap cMap) {
		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {
			@Override
			public ArrayList<Integer[]> call() throws Exception {
				Comparator<Integer[]> edge_cmp = generateContactMapEdgeComparator();
				ArrayList<Integer[]> edges = new ArrayList<>(cMap.edgeSet().size());
				for (Integer[] ent : cMap.edgeSet()) {
					int startTime = Population_Bridging.CONTACT_MAP_EDGE_START_TIME;
					int duration = Population_Bridging.CONTACT_MAP_EDGE_DURATION;
					while (duration < ent.length) {
						Integer[] e = new Integer[] { ent[Population_Bridging.CONTACT_MAP_EDGE_P1],
								ent[Population_Bridging.CONTACT_MAP_EDGE_P2], ent[startTime], ent[duration] };
						int key = Collections.binarySearch(edges, e, edge_cmp);
						
						if(key >= 0) {
							System.err.printf("Warning: Edge %s from %s already in list. Edge skipped.\n", 
									Arrays.toString(e), Arrays.toString(ent));
							
						}else {
							edges.add(~key, e);
						}
						startTime += 2;
						duration += 2;
					}
				}
				return edges;
			}
		};
		return callable;
	}

	public static final Callable<ArrayList<Integer[]>> generateContactEdgeArray(File cMap_file) {
		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {
			@Override
			public ArrayList<Integer[]> call() throws Exception {
				StringWriter cMap_str = new StringWriter();
				PrintWriter pWri = new PrintWriter(cMap_str);
				BufferedReader reader = new BufferedReader(new FileReader(cMap_file));
				String line;
				int minNumLine = 0;
				while ((line = reader.readLine()) != null) {
					pWri.println(line);
					minNumLine++;
				}
				pWri.close();
				reader.close();
				return generateContactEdgeArray(cMap_str.toString(), minNumLine).call();
			}
		};
		return callable;
	}

	public static final Callable<ArrayList<Integer[]>> generateContactEdgeArray(String cMap_str, int minSize) {
		Callable<ArrayList<Integer[]>> callable = new Callable<ArrayList<Integer[]>>() {

			@Override
			public ArrayList<Integer[]> call() throws Exception {
				ArrayList<Integer[]> edges = new ArrayList<>(minSize);
				BufferedReader lines = new BufferedReader(new StringReader(cMap_str));
				String line;
				Comparator<Integer[]> edge_cmp = generateContactMapEdgeComparator();
				
				while ((line = lines.readLine()) != null) {
					String[] ent = line.split(",");
					Integer p1 = Integer.parseInt(ent[Population_Bridging.CONTACT_MAP_EDGE_P1]);
					Integer p2 = Integer.parseInt(ent[Population_Bridging.CONTACT_MAP_EDGE_P2]);
					int startTime = Population_Bridging.CONTACT_MAP_EDGE_START_TIME;
					int duration = Population_Bridging.CONTACT_MAP_EDGE_DURATION;

					while (duration < ent.length) {
						Integer[] e = new Integer[] { p1, p2, Integer.parseInt(ent[startTime]),
								Integer.parseInt(ent[duration]), };

						int key = Collections.binarySearch(edges, e, edge_cmp);
						edges.add(~key, e);

						startTime += 2;
						duration += 2;
					}

				}

				return edges;
			}

		};

		return callable;
	}

	private static Comparator<Integer[]> generateContactMapEdgeComparator() {
		Comparator<Integer[]> edge_cmp = new Comparator<Integer[]>() {
			final int[] CMP_INDEX_ORDER = new int[] { Population_Bridging.CONTACT_MAP_EDGE_START_TIME,
					Population_Bridging.CONTACT_MAP_EDGE_DURATION, Population_Bridging.CONTACT_MAP_EDGE_P1,
					Population_Bridging.CONTACT_MAP_EDGE_P2, };

			@Override
			public int compare(Integer[] o1, Integer[] o2) {
				int cmp = 0;
				for (int cmpIndex : CMP_INDEX_ORDER) {
					cmp = Integer.compare(o1[cmpIndex], o2[cmpIndex]);
					if (cmp != 0) {
						break;
					}
				}
				return cmp;
			}
		};
		return edge_cmp;
	}

}
