package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel;
import sim.Simulation_ClusterModelGeneration;

public class Util_Combine_ContactMap {

	File dirSrcBase, dirSrcNew, dirRes;

	final Pattern PATTERN_CONTACT_MAP = Pattern.compile(
			Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

	public Util_Combine_ContactMap(File dirSrcBase, File dirSrcNew, File dirRes) {
		super();
		this.dirSrcBase = dirSrcBase;
		this.dirSrcNew = dirSrcNew;
		this.dirRes = dirRes;
	}

	public void combineMaps() throws FileNotFoundException, IOException {
		File[] base_map_arr = dirSrcBase.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return PATTERN_CONTACT_MAP.matcher(pathname.getName()).matches();
			}
		});

		for (File base_map : base_map_arr) {

			Matcher m = PATTERN_CONTACT_MAP.matcher(base_map.getName());
			m.matches();
			long mapSeed = Long.parseLong(m.group(1));

			File[] new_map_arr = dirSrcNew.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return Pattern.matches(Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP
							.replaceFirst("%d", Long.toString(mapSeed)).replaceFirst("%d", "(-{0,1}(?!0)\\\\d+)"),
							pathname.getName());
				}
			});

			if (new_map_arr.length > 0) {

				System.out.printf("Matching %d file(s) with contact map located in %s.\n", new_map_arr.length,
						base_map.getName());

				BufferedReader reader;
				StringWriter lines;
				PrintWriter pWri;
				String line;

				int maxEdgeStart = -1;

				reader = new BufferedReader(new FileReader(base_map));
				lines = new StringWriter();
				pWri = new PrintWriter(lines);

				while ((line = reader.readLine()) != null) {
					String[] ent = line.split(",");
					maxEdgeStart = Math.max(maxEdgeStart,
							Integer.parseInt(ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]));
					pWri.println(line);
				}
				reader.close();
				pWri.close();

				System.out.printf("Last edge in base map starts at %d.\n", maxEdgeStart);

				ContactMap base = ContactMap.ContactMapFromFullString(lines.toString());

				lines.close();

				for (File new_map : new_map_arr) {
					reader = new BufferedReader(new FileReader(new_map));
					while ((line = reader.readLine()) != null) {
						String[] ent = line.split(",");
						Integer[] newEdge = new Integer[ent.length];

						for (int i = 0; i < newEdge.length; i++) {
							newEdge[i] = Integer.parseInt(ent[i]);
						}

						Integer[] targetEdge = null;
						Set<Integer[]> edgeSet = base.getAllEdges(
								newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
								newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);

						if (edgeSet != null) {

							for (Integer[] testEdge : edgeSet) {
								if (testEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
										.equals(newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME])) {
									targetEdge = testEdge;
								}
							}
						}

						if (targetEdge != null) {
							targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = Math.max(
									targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION],
									newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]);
						} else {
							maxEdgeStart = Math.max(maxEdgeStart,
									newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
							targetEdge = newEdge;
							
							if(!base.containsVertex(targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1])) {
								base.addVertex(newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
							}
							
							if(!base.containsVertex(targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2])) {
								base.addVertex(newEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
							}
							
							
							base.addEdge(targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
									targetEdge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], targetEdge);

						}

					}
				}

				File targetFile = new File(dirRes, Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP
						.replaceFirst("%d", Long.toString(mapSeed)).replaceFirst("%d", Integer.toString(maxEdgeStart)));

				pWri = new PrintWriter(targetFile);
				pWri.println(base.toFullString());
				pWri.close();

				System.out.printf("Combine contact map %s generated.\n", targetFile.getName());

			}

		}

	}

}
