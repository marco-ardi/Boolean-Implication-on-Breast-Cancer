import java.util.Vector;

public class StepMiner
{
    public static void main(String[] args) throws Exception
    {
        String inputFile=null;
        double delta=0.5;
        String outputFile="expr.txt";

        //Read input parameters
        int i=0;
        for (i=0;i<args.length;i++)
        {
            switch (args[i])
            {
                case "-i" -> inputFile = args[++i];
                case "-d" -> delta = Double.parseDouble(args[++i]);
                case "-o" -> outputFile = args[++i];
                default -> {
                    System.out.println("Error! Unrecognizable command '" + args[i]);
                    printHelp();
                    System.exit(1);
                }
            }
        }

        //Error in case network file is missing
        if(inputFile==null)
        {
            System.out.println("Error! No expression matrix has been specified!\n");
            printHelp();
            System.exit(1);
        }

        System.out.println("Discretizing expression matrix...");

        //Read expression data
        FileManager fm=new FileManager();
        Vector<double[]> mapExpressions=fm.readExpressionFile(inputFile);

        //Start boolean analysis
        BooleanAnalyzer sa=new BooleanAnalyzer();

        //Run StepMiner to discretize expression values
        Vector<int[]> listDiscretizedValues=sa.discretizeMatrix(mapExpressions,delta);

        //Write output results to file
        fm.writeDiscretizedMatrix(listDiscretizedValues,outputFile);

        System.out.println("Done!");
    }

    private static void printHelp()
    {
        String help = "Usage: java StepMiner -m <matrixFile> [-d <delta> ";
        help+="-o <outputFile>]\n\n";
        help+="REQUIRED PARAMETERS:\n";
        help+="-i\tExpression matrix file\n\n";
        help+="OPTIONAL PARAMETERS:\n";
        help+="-d\tDelta threshold around center for low and high values (default=0.5)\n";
        help+="-o\tOutput file where discretized expression matrix will be saved (default=expr.txt)\n\n";
        System.out.println(help);
    }
}
