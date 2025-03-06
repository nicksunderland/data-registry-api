import json
import subprocess
import boto3
import click


def download_file_from_s3(s3_path):
    s3 = boto3.client('s3')
    bucket, key = s3_path.replace("s3://", "").split("/", 1)
    remote_file_name = key.split('/')[-1]
    s3.download_file(bucket, key, remote_file_name)
    return remote_file_name


def upload_file_to_s3(file_name, file_guid, extra_args={'ContentType': 'image/png'}):
    s3 = boto3.client('s3')
    bucket = "hermes-qc"
    key = f"images/{file_guid}/" + file_name
    s3.upload_file(file_name, bucket, key, ExtraArgs=extra_args)


