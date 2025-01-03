import os
import json
import requests
import multiprocess
import zipfile

def get_data_from_dataverse(data):
    file_url = data.get("url")
    name_file = data.get("name_file")
    resp = requests.get(file_url)
    with open(name_file, "wb") as f:
        f.write(resp.content)
    
    # Decompress the file
    with zipfile.ZipFile(name_file, 'r') as zip_ref:
        zip_ref.extractall("files/xtc_files")

def download_xtc_files():
    folder_already_exist=True
    if not os.path.exists("files/xtc_files"):
        os.mkdir("files/xtc_files")
        folder_already_exist=False

    if not folder_already_exist:
        with open("./files/json/dataverse.json", "r") as file:
            urls = json.load(file)

        print("Your download is starting...")

        pool = multiprocess.Pool(processes=len(urls))
        pool.map(get_data_from_dataverse, urls)

        print("You have successfully downloaded and decompressed the data from Dataverse.")
    else:
        print("You have download the data previously")

if __name__ == "__main__":
    download_xtc_files()
