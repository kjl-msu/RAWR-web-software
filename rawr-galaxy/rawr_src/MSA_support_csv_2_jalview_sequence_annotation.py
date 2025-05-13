import argparse
import pandas as pd
def support_csv_2_jalview(MSA_support_csv,color):
    cols=["columnIndex","supportValue"]
    try:
        df = pd.read_csv(MSA_support_csv, usecols=cols, sep=",")
    except ValueError:
        print("File does not have the right column names")
    df = df.groupby(by="columnIndex").mean()
    df = df.reset_index()
    numcols = int(df.iloc[-1]["columnIndex"])
    jalview_bargraph_str="JALVIEW_ANNOTATION\nBAR_GRAPH\tSupport\t"
    for i in range(0, numcols):
        try:
            support = str(df.iloc[i]["supportValue"])
        except IndexError:
            support = str("0")
        jalview_bargraph_str = jalview_bargraph_str + support+","+support+"|"
    jalview_bargraph_str = jalview_bargraph_str + "\nCOLOUR\tSupport\t"+color
    with open(MSA_support_csv+"_jalview_annotation.txt","w") as out:
        out.write(jalview_bargraph_str)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert RAWR software MSA support CSV file to JalView sequence features annotation file, where the support values are displayed as bar graphs.')
    parser.add_argument('MSA_support_csv', type=str, help='<CSV> File output from running MSA support with RAWR software.')
    parser.add_argument('-color', type=str, help='<optional> Set JalView bar graph color', default="pink")

    args = parser.parse_args()
    support_csv_2_jalview(args.MSA_support_csv,args.color)
