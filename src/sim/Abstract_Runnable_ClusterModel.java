package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.Callable;

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
	
	

}
