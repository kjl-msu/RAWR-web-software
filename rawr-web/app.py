# Available under GPLv3 license
import os
import time
from flask import Flask, request, redirect, url_for, render_template, make_response, send_file, send_from_directory
from werkzeug.utils import secure_filename
from src import sampleSeq, seqs
import subprocess
from flask_mail import Mail, Message
import sys

app = Flask(__name__, static_url_path='')

UPLOAD_FOLDER = 'upload'
RESULT_FOLDER = 'result'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULT_FOLDER'] = RESULT_FOLDER

"""
Flask with gmail as server sender. (You must provide a mail server to send emails, and the easiest way is to use an existing service.)

For Gmail, please enable 2-step verification and then generate an app password to be used here as MAIL_PASSWORD. Example:
        MAIL_SERVER = 'smtp.gmail.com'
        MAIL_PORT = 465
        MAIL_USE_SSL = True
        MAIL_USERNAME = 'username@gmail.com'
        MAIL_PASSWORD = 'app password generated by Gmail'
"""
app.config.update(
    MAIL_SERVER='smtp.gmail.com',
    MAIL_PORT=465,
    MAIL_USE_SSL=True,
    MAIL_USERNAME='username@gmail.com',
    MAIL_PASSWORD='app password generated by Gmail'
)
mail = Mail(app)
basedir = os.path.abspath(os.path.dirname(__file__))
ALLOWED_EXTENSIONS = {'fasta','fas','fa'} #set(['fasta'])
# scriptDir to where the server directory should be, must have forward slash at the end 
scriptDir = os.path.dirname(os.path.abspath(__file__))+"/src/"

def allowed_file(filename):
    """
    :param filename: file name.
    :return: return if the file extension is valid for upload.
    """
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


def runShellCmd(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate()[0]
    return output


def writeSampleSeqAndIndex(sampleSeqData, sampleIndex, outputPath, n):
    """
    Write sampled sequences and index to output file.
    :param sampleSeqData: pandas.DataFrame, sequence data.
    :param sampleIndex: numpy.Array, index data.
    :param outputPath: string, output prefix.
    :param n: int, sample number.
    :return:
    """
    with open(os.path.join(outputPath,str(n) + ".seq.fasta"), 'w') as outf:
        for i in sampleSeqData.index:
            outf.write('>' + sampleSeqData.loc[i].name + '\n')
            outf.write(''.join(sampleSeqData.loc[i]).replace('-', '') +
                       '\n')
    with open(os.path.join(outputPath, str(n) + ".index"), "w") as outf:
        for idx in sampleIndex:
            outf.write(str(idx) + "\n")


def generateSampleSeq(alndata, outputPath, samplenum,
                      sampleMethod, methodParameters):
    """
    Generate SERES or RAWR resampled sequences.
    :param alndata: pandas.DataFrame, input alignment data.
    :param outputPath: string, output prefix.
    :param samplenum: int, number of sampled replicates.
    :param sampleMethod: string, "SERES" or "RAWR".
    :param methodParameters: list, parameters for sampling algorithm.
    :return: None.
    """
    if sampleMethod == "RAWR":
        reverseRate = methodParameters[0]
        for n in range(1, samplenum + 1):
            sampleIndex, sampleSeqData = sampleSeq.rawrSample(
                alndata, float(reverseRate))
            writeSampleSeqAndIndex(sampleSeqData, sampleIndex, outputPath, n)

    elif sampleMethod == "SERES":
        anchorLen, anchorNum, reverseRate = methodParameters
        anchorLen = "default" if not anchorLen.isdigit() else int(anchorLen)
        anchorNum = "default" if not anchorNum.isdigit() else int(anchorNum)
        barrier = sampleSeq.getAnchor(
            alndata, int(anchorLen), int(anchorNum))
        for n in range(1, samplenum + 1):
            sampleIndex, sampleSeqData = sampleSeq.seresSample(
                alndata, int(anchorLen), int(anchorNum), float(reverseRate), barrier)
            writeSampleSeqAndIndex(sampleSeqData, sampleIndex, outputPath, n)
    else:
        print(
            'Please assign a sampling method: RAWR, SERES'
        )


def sampleSequence(workDir, sampleNum, sampleMethod,
                   methodParameters, supportType, packupParams):
    alnfile = os.path.join(workDir, "alignment.fasta")
    outPrefix = os.path.join(workDir, "samples")
    alndata = seqs.getAlnData(alnfile)
    sampleNum = int(sampleNum)
    generateSampleSeq(
        alndata,
        outPrefix,
        sampleNum,
        sampleMethod,
        methodParameters)
    if supportType == "MSA":
        cmd = " ".join(["bash", scriptDir + "calMSASupport.sh",
            workDir, str(sampleNum), scriptDir])
        print(cmd)
        runShellCmd(cmd)
    elif supportType == "Tree":
        cmd = " ".join(["bash", scriptDir + "calTreeSupport.sh",
            workDir, str(sampleNum), scriptDir])
        print(cmd)
        runShellCmd(cmd)
    cmd = " ".join(["bash", scriptDir + "packupForDownload.sh"] + packupParams)
    print(cmd)
    runShellCmd(cmd)

def mail_results(filename, emailAddress, supportFile, jalviewFile, sampleFile, TreeSupportFigure=None):
    ### email user asynchronously with link to download these files
    msg_title = 'RAWR resampling results '+filename
    msg_sender = 'Sender '+app.config['MAIL_USERNAME']
    msg_recipients = [emailAddress]
    msg_body = "Your RAWR v1.0 job has finished. \nPlease find attached a Support File, a Sampling File"
    mail_files = [supportFile, sampleFile]
    if TreeSupportFigure is not None:
        mail_files.append(TreeSupportFigure)
        msg_body += " and a Tree Support Figure.\n"
    else:
        mail_files.append(jalviewFile)
        msg_body += " and a JalView annotation file. \n"
    msg = Message(msg_title,
                  sender=msg_sender,
                  recipients=msg_recipients,
                  body=msg_body)
    for file in mail_files:
        print(file)
        with app.open_resource(file) as fp:
            msg.attach(file,"text/plain", fp.read())
    print(msg.body)

    print("sending mail")
    mail.send(msg)
    print("mail sent")

@app.route('/', methods=['GET', 'POST'])
def main_page():
    if request.method == 'POST':
        uploadDir = os.path.join(
            basedir, app.config['UPLOAD_FOLDER'])
        resultDir = os.path.join(
            basedir, app.config['RESULT_FOLDER'])
        if not os.path.exists(uploadDir):
            os.makedirs(uploadDir)
        if not os.path.exists(resultDir):
            os.makedirs(resultDir)
        file = request.files['alignmentFile']
        filename = secure_filename(file.filename)
        if file and allowed_file(filename):
            unix_time = str(int(time.time()))
            workDir = os.path.join(uploadDir, unix_time)
            os.makedirs(workDir)
            sampleDir = os.path.join(workDir, "samples")
            os.makedirs(sampleDir)
            alnfile = os.path.join(workDir, "alignment.fasta")
            file.save(alnfile)

            reverseRate = request.form["reverseRate"]
            sampleNum = request.form["sampleNum"]
            sampleMethod = request.form["algorithm"]
            supportType = request.form["supportType"]
            emailAddress = request.form["emailAddress"]

            if supportType == "Tree":
                phyloTree = request.form["phyloTree"]
                with open(os.path.join(workDir, "input.tree"), "w") as outf:
                    outf.write(phyloTree + "\n")
                supportFile = "tree.support.txt"
                TreeSupportFigure = os.path.join(workDir, "tree.support.png")
            else:
                supportFile = "MSA_Support.csv"
                TreeSupportFigure = None
            app.config['SUPPORT_FILE'] = supportFile
            packupParams = [basedir, unix_time, supportFile]
            if sampleMethod == "RAWR":
                sampleSequence(workDir, sampleNum, sampleMethod,
                               [reverseRate], supportType, packupParams)
            elif sampleMethod == "SERES":
                anchorLen = request.form["anchorLen"]
                anchorNum = request.form["anchorNum"]
                sampleSequence(workDir, sampleNum, sampleMethod,
                               [anchorLen, anchorNum, reverseRate], supportType, packupParams)

            supportFile_result = os.path.join('upload',str(unix_time), app.config['SUPPORT_FILE'])
            jalviewFile_result = os.path.join('upload',str(unix_time), app.config['SUPPORT_FILE'] +"_jalview_annotation.txt")
            sampleFile_result = os.path.join('upload',str(unix_time), "samples.tar.gz")
            mail_results(unix_time, emailAddress, supportFile_result,jalviewFile_result, sampleFile_result, TreeSupportFigure)
            return redirect(url_for("result", filename=unix_time))
    return render_template('index.html')


@app.route('/result/<filename>')     # server's result folder. will autoreload server if script is changed
def result(filename):
    supportFile = filename + "." + app.config['SUPPORT_FILE']
    jalviewFile = filename + "." + app.config['SUPPORT_FILE'] + "_jalview_annotation.txt"
    sampleFile = filename + "." + "samples.tar.gz"
    downloads = {"supportFile": supportFile, "sampleFile": sampleFile}

    if app.config['SUPPORT_FILE'] == "tree.support.txt":
        downloads["figFile"] = filename + ".tree.support.png"
        print(downloads)
        return render_template('treeResult.html', file=downloads)
    elif app.config['SUPPORT_FILE'] == "MSA_Support.csv":
        downloads["jalviewFile"] = jalviewFile
        return render_template('msaResult.html', file=downloads)

@app.route('/results_img/<filename>')
def results_img(filename):
    resultDir = os.path.join(basedir, app.config['RESULT_FOLDER'])
    return send_from_directory(resultDir, filename)

@app.route("/download/<filename>", methods=['GET'])
def download(filename):
    resultDir = os.path.join(
        basedir, app.config['RESULT_FOLDER'])
    response = make_response(send_from_directory(
        resultDir, filename.encode('utf-8').decode('utf-8'), as_attachment=True))
    response.headers["Content-Disposition"] = "attachment; filename={}".format(
        filename.encode().decode('latin-1'))

    return response


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
