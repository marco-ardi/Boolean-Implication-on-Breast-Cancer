import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

public class BooleanAnalyzer
{
    public Vector<int[]> discretizeMatrix(Vector<double[]> listValues, double gap)
    {
        Vector<int[]> listDiscretizedValues=new Vector<>();
        for(double[] originalValues : listValues)
        {
            double[] values=Arrays.copyOf(originalValues,originalValues.length);
            Arrays.sort(values);
            int bestPos=-1;
            double minSSE=Double.MAX_VALUE;
            for(int i=0;i<values.length;i++)
            {
                double currentSSE=0.0;
                double leftSum=0.0;
                for(int j=0;j<=i;j++)
                    leftSum+=values[j];
                leftSum/=(i+1);
                for(int j=0;j<=i;j++)
                    currentSSE+=Math.pow(values[i]-leftSum,2);
                double rightSum=0.0;
                for(int j=i+1;j<values.length;j++)
                    rightSum+=values[j];
                rightSum/=(values.length-i-1);
                for(int j=i+1;j<values.length;j++)
                    currentSSE+=Math.pow(values[i]-rightSum,2);
                //System.out.println(currentSSE);
                if(currentSSE<minSSE)
                {
                    minSSE=currentSSE;
                    bestPos=i;
                }
            }
            //System.out.println(bestPos);
            double threshold=values[bestPos];
            //System.out.println(threshold);
            double lowerBound=threshold-gap;
            double upperBound=threshold+gap;
            int[] discretizedValues=new int[originalValues.length];
            for(int i=0;i<originalValues.length;i++)
            {
                if(originalValues[i]>upperBound)
                    discretizedValues[i]=1;
                else if(originalValues[i]<lowerBound)
                    discretizedValues[i]=-1;
            }
            listDiscretizedValues.add(discretizedValues);
        }
        return listDiscretizedValues;
    }

    private int[][] getQuadrantCounts(String geneA, String geneB, HashMap<String,int[]> mapValues)
    {
        int[][] quadrantCounts=new int[2][2];
        int[] valuesA=mapValues.get(geneA);
        int[] valuesB=mapValues.get(geneB);
        for(int i=0;i<valuesA.length;i++)
        {
            if(valuesA[i]==-1 && valuesB[i]==-1)
                quadrantCounts[0][0]++;
            else if(valuesA[i]==-1 && valuesB[i]==1)
                quadrantCounts[0][1]++;
            else if(valuesA[i]==1 && valuesB[i]==-1)
                quadrantCounts[1][0]++;
            else if(valuesA[i]==1 && valuesB[i]==1)
                quadrantCounts[1][1]++;
        }
        return quadrantCounts;
    }

    public Vector<Vector<String>> getImplications(String geneA, String type, HashMap<String,int[]> mapValues,
                                                  double statThresh, double pvalThresh)
    {
        Vector<Vector<String>> listImplications=new Vector<>();
        for(String geneB : mapValues.keySet())
        {
            if(!geneB.equals(geneA))
                listImplications.addAll(getImplications(geneA,geneB,type,mapValues,statThresh,pvalThresh));
        }
        return listImplications;
    }

    public Vector<Vector<String>> getImplications(String geneA, String geneB, String type, HashMap<String,int[]> mapValues,
                                                  double statThresh, double pvalThresh)
    {
        Vector<Vector<String>> listImplications=new Vector<>();
        int[][] quadrantCounts=getQuadrantCounts(geneA,geneB,mapValues);
        if(type.equals("any"))
        {
            listImplications.addAll(getImplications(geneA,geneB,"low-low",quadrantCounts,statThresh,pvalThresh));
            listImplications.addAll(getImplications(geneA,geneB,"low-high",quadrantCounts,statThresh,pvalThresh));
            listImplications.addAll(getImplications(geneA,geneB,"high-low",quadrantCounts,statThresh,pvalThresh));
            listImplications.addAll(getImplications(geneA,geneB,"high-high",quadrantCounts,statThresh,pvalThresh));
            listImplications.addAll(getImplications(geneA,geneB,"equivalent",quadrantCounts,statThresh,pvalThresh));
            listImplications.addAll(getImplications(geneA,geneB,"opposite",quadrantCounts,statThresh,pvalThresh));
        }
        else
            listImplications.addAll(getImplications(geneA,geneB,type,quadrantCounts,statThresh,pvalThresh));
        return listImplications;
    }

    private Vector<Vector<String>> getImplications(String geneA, String geneB, String type, int[][] quadrantCounts,
                                                  double statThresh, double pvalThresh)
    {
        Vector<Vector<String>> listImplications=new Vector<>();
        int nTotal=quadrantCounts[0][0]+quadrantCounts[0][1]+quadrantCounts[1][0]+quadrantCounts[1][1];
        int nFirstLow=quadrantCounts[0][0]+quadrantCounts[0][1];
        int nFirstHigh=quadrantCounts[1][0]+quadrantCounts[1][1];
        int nSecondHigh=quadrantCounts[0][1]+quadrantCounts[1][1];
        int nSecondLow=quadrantCounts[0][0]+quadrantCounts[1][0];
        if(nTotal>0)
        {
            if(type.equals("low-low"))
            {
                if(nFirstLow>0 && nSecondHigh>0)
                {
                    double nExpected=((double)nFirstLow*nSecondHigh)/nTotal;
                    double statistic=(nExpected-quadrantCounts[0][1])/Math.sqrt(nExpected);
                    double pval=((((double)quadrantCounts[0][1])/nFirstLow)+(((double)quadrantCounts[0][1])/nSecondHigh))/2;
                    if(statistic>=statThresh && pval<=pvalThresh)
                    {
                        Vector<String> newImplication=new Vector<>();
                        newImplication.add(geneA);
                        newImplication.add(geneB);
                        newImplication.add(type);
                        newImplication.add(String.valueOf(statistic));
                        newImplication.add(String.valueOf(pval));
                        listImplications.add(newImplication);
                    }
                }
            }
            else if(type.equals("low-high"))
            {
                if(nFirstLow>0 && nSecondLow>0)
                {
                    double nExpected=((double)nFirstLow*nSecondLow)/nTotal;
                    double statistic=(nExpected-quadrantCounts[0][0])/Math.sqrt(nExpected);
                    double pval=((((double)quadrantCounts[0][0])/nFirstLow)+(((double)quadrantCounts[0][0])/nSecondLow))/2;
                    if(statistic>=statThresh && pval<=pvalThresh)
                    {
                        Vector<String> newImplication=new Vector<>();
                        newImplication.add(geneA);
                        newImplication.add(geneB);
                        newImplication.add(type);
                        newImplication.add(String.valueOf(statistic));
                        newImplication.add(String.valueOf(pval));
                        listImplications.add(newImplication);
                    }
                }
            }
            else if(type.equals("high-low"))
            {
                if(nFirstHigh>0 && nSecondHigh>0)
                {
                    double nExpected=((double)nFirstHigh*nSecondHigh)/nTotal;
                    double statistic=(nExpected-quadrantCounts[1][1])/Math.sqrt(nExpected);
                    double pval=((((double)quadrantCounts[1][1])/nFirstHigh)+(((double)quadrantCounts[1][1])/nSecondHigh))/2;
                    if(statistic>=statThresh && pval<=pvalThresh)
                    {
                        Vector<String> newImplication=new Vector<>();
                        newImplication.add(geneA);
                        newImplication.add(geneB);
                        newImplication.add(type);
                        newImplication.add(String.valueOf(statistic));
                        newImplication.add(String.valueOf(pval));
                        listImplications.add(newImplication);
                    }
                }
            }
            else if(type.equals("high-high"))
            {
                if(nFirstHigh>0 && nSecondLow>0)
                {
                    double nExpected=((double)nFirstHigh*nSecondLow)/nTotal;
                    double statistic=(nExpected-quadrantCounts[1][0])/Math.sqrt(nExpected);
                    double pval=((((double)quadrantCounts[1][0])/nFirstHigh)+(((double)quadrantCounts[1][0])/nSecondLow))/2;
                    if(statistic>=statThresh && pval<=pvalThresh)
                    {
                        Vector<String> newImplication=new Vector<>();
                        newImplication.add(geneA);
                        newImplication.add(geneB);
                        newImplication.add(type);
                        newImplication.add(String.valueOf(statistic));
                        newImplication.add(String.valueOf(pval));
                        listImplications.add(newImplication);
                    }
                }
            }
            else if(type.equals("equivalent"))
            {
                if(nFirstLow>0 && nSecondHigh>0)
                {
                    double nExpected1=((double)nFirstLow*nSecondHigh)/nTotal;
                    double statistic1=(nExpected1-quadrantCounts[0][1])/Math.sqrt(nExpected1);
                    double pval1=((((double)quadrantCounts[0][1])/nFirstLow)+(((double)quadrantCounts[0][1])/nSecondHigh))/2;
                    if(nFirstHigh>0 && nSecondLow>0)
                    {
                        double nExpected2=((double)nFirstHigh*nSecondLow)/nTotal;
                        double statistic2=(nExpected2-quadrantCounts[1][0])/Math.sqrt(nExpected2);
                        double pval2=((((double)quadrantCounts[1][0])/nFirstHigh)+(((double)quadrantCounts[1][0])/nSecondLow))/2;
                        if(statistic1>=statThresh && pval1<=pvalThresh && statistic2>=statThresh && pval2<=pvalThresh)
                        {
                            Vector<String> newImplication=new Vector<>();
                            newImplication.add(geneA);
                            newImplication.add(geneB);
                            newImplication.add(type);
                            newImplication.add(String.valueOf(statistic1));
                            newImplication.add(String.valueOf(pval1));
                            newImplication.add(String.valueOf(statistic2));
                            newImplication.add(String.valueOf(pval2));
                            listImplications.add(newImplication);
                        }
                    }
                }
            }
            else
            {
                if(nFirstLow>0 && nSecondLow>0)
                {
                    double nExpected1=((double)nFirstLow*nSecondLow)/nTotal;
                    double statistic1=(nExpected1-quadrantCounts[0][0])/Math.sqrt(nExpected1);
                    double pval1=((((double)quadrantCounts[0][0])/nFirstLow)+(((double)quadrantCounts[0][0])/nSecondLow))/2;
                    if(nFirstHigh>0 && nSecondHigh>0)
                    {
                        double nExpected2=((double)nFirstHigh*nSecondHigh)/nTotal;
                        double statistic2=(nExpected2-quadrantCounts[1][1])/Math.sqrt(nExpected2);
                        double pval2=((((double)quadrantCounts[1][1])/nFirstHigh)+(((double)quadrantCounts[1][1])/nSecondHigh))/2;
                        if(statistic1>=statThresh && pval1<=pvalThresh && statistic2>=statThresh && pval2<=pvalThresh)
                        {
                            Vector<String> newImplication=new Vector<>();
                            newImplication.add(geneA);
                            newImplication.add(geneB);
                            newImplication.add(type);
                            newImplication.add(String.valueOf(statistic1));
                            newImplication.add(String.valueOf(pval1));
                            newImplication.add(String.valueOf(statistic2));
                            newImplication.add(String.valueOf(pval2));
                            listImplications.add(newImplication);
                        }
                    }
                }
            }
        }
        return listImplications;
    }
}
