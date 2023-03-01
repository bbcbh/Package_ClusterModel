package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import population.Population_Bridging;
import sim.Simulation_ClusterModelGeneration;

public class Util_Analyse_ContactMap_Outputs {

	File basedir;
	int numCat;
	int keepVal;

	final Pattern REGEX_FIRST_LINE = Pattern.compile("Seed = (-?\\d+)");
	final Pattern REGEX_OUTPUT_FILENAME = Pattern
			.compile(Simulation_ClusterModelGeneration.FILENAME_FORMAT_OUTPUT.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));
	final Pattern REGEX_INDEX_LINE = Pattern
			.compile("# partners in last 12 months -  Day %d".replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));
	final Pattern REGEX_ENTRY_LINE;

	public Util_Analyse_ContactMap_Outputs(String baseDirPath, int numCat, int keepVal) {
		basedir = new File(baseDirPath);
		this.numCat = numCat;
		this.keepVal = keepVal;
		StringBuilder entryStr = new StringBuilder();
		entryStr.append(" (\\d): \\[(\\d+)");
		for (int i = 1; i < numCat; i++) {
			entryStr.append(", (\\d+)");
		}
		entryStr.append("\\]");
		REGEX_ENTRY_LINE = Pattern.compile(entryStr.toString());
	}

	public void analyse_output() throws FileNotFoundException, IOException {
		File[] outputFiles = basedir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return REGEX_OUTPUT_FILENAME.matcher(pathname.getName()).matches();
			}
		});

		System.out.printf("%d output file(s) found in %s.\n", outputFiles.length, basedir.getAbsolutePath());

		@SuppressWarnings("unchecked")
		// Key = day, value=double[] {% in category for sim #1}
		HashMap<Integer, ArrayList<Double>>[][] entry_collection = new HashMap[Population_Bridging.LENGTH_GENDER][numCat];

		for (File outputFile : outputFiles) {

			BufferedReader reader = new BufferedReader(new FileReader(outputFile));

			ArrayList<String> linesCollection = new ArrayList<>();
			String line;
			int lineNum = 1;

			while ((line = reader.readLine()) != null) {
				if (REGEX_FIRST_LINE.matcher(line).matches()) {
					// Reset
					linesCollection.clear();
					if (lineNum != 1) {
						System.out.printf("  Warning: Lines reading restart to line %d for %s\n", lineNum,
								outputFile.getName());
						break;
					}

				}
				if (line.length() > 0) {
					linesCollection.add(line);
				}
				lineNum++;
			}
			reader.close();

			Iterator<String> lines = linesCollection.iterator();
			Matcher m;

			while (lines.hasNext()) {
				line = lines.next();
				m = REGEX_INDEX_LINE.matcher(line);
				if (m.matches()) {
					int day = Integer.parseInt(m.group(1));

					int[][] oneDayEntry = new int[entry_collection.length][numCat];

					// Read and check entries
					try {
						for (int g = 0; g < entry_collection.length; g++) {
							line = lines.next();
							m = REGEX_ENTRY_LINE.matcher(line);
							if (m.matches()) {
								int gender = Integer.parseInt(m.group(1));
								if (gender != g) {
									System.err.printf(
											"Warning: Gender line mismatches for line \"%s\" with gender #%d in %s\n",
											line, g, outputFile.getName());
								}

								for (int c = 0; c < oneDayEntry[gender].length; c++) {
									oneDayEntry[gender][c] = Integer.parseInt(m.group(c + 2));
								}
							} else {
								oneDayEntry = null;
								break;
							}
						}

					} catch (Exception ex) {
						ex.printStackTrace(System.err);
						oneDayEntry = null;
					}
					if (oneDayEntry != null) {
						for (int g = 0; g < oneDayEntry.length; g++) {
							double total = 0;
							for (int c = 0; c < oneDayEntry[g].length; c++) {
								total += oneDayEntry[g][c];
							}

							if (total != 0) {
								for (int c = 0; c < oneDayEntry[g].length; c++) {
									double proportion = oneDayEntry[g][c] / total;

									HashMap<Integer, ArrayList<Double>> map = entry_collection[g][c];
									if (map == null) {
										map = new HashMap<>();
										entry_collection[g][c] = map;
									}
									ArrayList<Double> ent = map.get(day);
									if (ent == null) {
										ent = new ArrayList<>();
										map.put(day, ent);
									}
									ent.add(proportion);
								}
							}
						}
					}
				}
			}
			System.out.printf("Analysis of %s completed.\n", outputFile.getName());
		}

		System.out.println();
		for (int g = 0; g < entry_collection.length; g++) {			
			for (int c = 0; c < entry_collection[g].length; c++) {			
				HashMap<Integer, ArrayList<Double>> map = entry_collection[g][c];
				if (map != null) {					
					System.out.printf("Gender %d, Cat %d:\n", g, c);
					Integer[] days = map.keySet().toArray(new Integer[map.size()]);
					DescriptiveStatistics statAllDays = new DescriptiveStatistics(days.length);

					Arrays.sort(days);
					for (int i = Math.max(0, days.length - keepVal); i < days.length; i++) {												
						Integer day = days[i];
						ArrayList<Double> vals = map.get(day);
						DescriptiveStatistics stat = new DescriptiveStatistics(vals.size());
						for (Double val : vals) {
							stat.addValue(val);
							statAllDays.addValue(val);
						}
						System.out.printf(" %7d: %.3f, (%.3f - %.3f)\n", day, stat.getPercentile(50), stat.getPercentile(25),
								stat.getPercentile(75));
					}

					System.out.printf(" Last %d: %.3f, (%.3f - %.3f)\n", Math.min(days.length, keepVal), statAllDays.getPercentile(50),
							statAllDays.getPercentile(25), statAllDays.getPercentile(75));
				}

			}
		}

	}

}
