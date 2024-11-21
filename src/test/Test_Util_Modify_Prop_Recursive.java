package test;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Properties;
import java.util.regex.Pattern;

public class Test_Util_Modify_Prop_Recursive {

	final static Pattern pattern_xml_filename = Pattern.compile("(.*)_simSpecificSim.prop");
	final static String[][] mod_entries = new String[][] {
			new String[] { "PROP_PEP_EFFICACY", "[-0.77,-0.22,-0.78]" }, };

	final static FileFilter filesearch_filter = new FileFilter() {
		@Override
		public boolean accept(File pathname) {
			return pathname.isDirectory() || pattern_xml_filename.matcher(pathname.getName()).find();
		}
	};

	public static void main(String[] args) {
		File baseDir = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Gen");
		long tic = System.currentTimeMillis();
		modify_prop(baseDir);
		System.out.printf("Generation of new files completed. Time req. = %.3fs",
				(System.currentTimeMillis() - tic) / 1000f);
	}

	protected static void modify_prop(File baseDir) {
		File[] fileList = baseDir.listFiles(filesearch_filter);
		for (File f : fileList) {
			if (f.isDirectory()) {
				modify_prop(f);
			} else {				
				Properties prop = new Properties();
				try {
					prop.loadFromXML(new FileInputStream(f));
					long timeStamp = System.currentTimeMillis();
					String backupfName = String.format("%d_%s", timeStamp, f.getName());

					String comment = String.format(
							"Generated at %1$tF %1$tT based on %2$s\n with modiflied entries %3$s", timeStamp,
							backupfName, Arrays.deepToString(mod_entries));

					for (String[] ent : mod_entries) {
						prop.put(ent[0], ent[1]);
					}

					// Make a backup copy
					Files.copy(f.toPath(), new File(baseDir, backupfName).toPath());

					prop.storeToXML(new FileOutputStream(f), comment);

				} catch (IOException e) {
					e.printStackTrace(System.err);
				}

			}
		}

	}

}
