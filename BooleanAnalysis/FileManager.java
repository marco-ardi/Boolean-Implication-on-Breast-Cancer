import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;

public class FileManager
{
    private String[] listSamples;
    private final Vector<String> listGenes;

    public FileManager()
    {
        listGenes=new Vector<>();
    }

    public Vector<double[]> readExpressionFile(String file) throws Exception
    {
        Vector<double[]> mapExpressions=new Vector<>();
        BufferedReader br=new BufferedReader(new FileReader(file));
        String str;
        listSamples=br.readLine().split("\t");
        while((str=br.readLine())!=null)
        {
            String[] split=str.split("\t");
            String gene=split[0];
            double[] values=new double[listSamples.length];
            for(int i=1;i<split.length;i++)
                values[i-1]=Double.parseDouble(split[i]);
            mapExpressions.add(values);
            listGenes.add(gene);
        }
        br.close();
        return mapExpressions;
    }

    public HashMap<String,int[]> readDiscretizedExpressionFile(String file) throws Exception
    {
        HashMap<String,int[]> mapExpressions=new HashMap<>();
        BufferedReader br=new BufferedReader(new FileReader(file));
        String str;
        listSamples=br.readLine().split("\t");
        while((str=br.readLine())!=null)
        {
            String[] split=str.split("\t");
            String gene=split[0];
            int[] values=new int[listSamples.length];
            for(int i=1;i<split.length;i++)
                values[i-1]=Integer.parseInt(split[i]);
            mapExpressions.put(gene,values);
        }
        br.close();
        return mapExpressions;
    }

    public void initImplicationFile(String file) throws Exception
    {
        BufferedWriter bw=new BufferedWriter(new FileWriter(file));
        bw.write("Implication\tStatistic(s)\tP-value(s)\n");
        bw.close();
    }

    public void writeImplications(Vector<Vector<String>> listImplications, String file) throws Exception
    {
        BufferedWriter bw=new BufferedWriter(new FileWriter(file,true));
        for (Vector<String> implication : listImplications)
        {
            String geneA = implication.get(0);
            String geneB = implication.get(1);
            String type = implication.get(2);
            if (type.equals("low-low"))
                bw.write(geneA + " low => " + geneB + " low\t"+implication.get(3)+"\t"+implication.get(4)+"\n");
            else if (type.equals("low-high"))
                bw.write(geneA + " low => " + geneB + " high\t"+implication.get(3)+"\t"+implication.get(4)+"\n");
            else if (type.equals("high-low"))
                bw.write(geneA + " high => " + geneB + " low\t"+implication.get(3)+"\t"+implication.get(4)+"\n");
            else if (type.equals("high-high"))
                bw.write(geneA + " high => " + geneB + " high\t"+implication.get(3)+"\t"+implication.get(4)+"\n");
            else if (type.equals("equivalent"))
            {
                bw.write(geneA + " low => " + geneB + " low AND " + geneA + " high => " + geneB + " high\t");
                bw.write(implication.get(3)+","+implication.get(5)+"\t"+implication.get(4)+","+implication.get(6)+"\n");
            }
            else
            {
                bw.write(geneA + " low => " + geneB + " high AND " + geneA + " high => " + geneB + " low\n");
                bw.write(implication.get(3)+","+implication.get(5)+"\t"+implication.get(4)+","+implication.get(6)+"\n");
            }
        }
        bw.close();
    }

    public void writeDiscretizedMatrix(Vector<int[]> listValues, String file) throws Exception
    {
        BufferedWriter bw=new BufferedWriter(new FileWriter(file));
        int i=0;
        for(i=0;i<listSamples.length-1;i++)
            bw.write(listSamples[i]+"\t");
        bw.write(listSamples[i]+"\n");
        for(i=0;i<listValues.size();i++)
        {
            bw.write(listGenes.get(i));
            int[] values=listValues.get(i);
            for(int value : values)
                bw.write("\t" + value);
            bw.write("\n");
        }
        bw.close();
    }

    public String[] getListSamples()
    {
        return listSamples;
    }

    public Vector<String> getListGenes()
    {
        return listGenes;
    }

}
