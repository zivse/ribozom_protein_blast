import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file_path', help='foo help')
parser.add_argument('--directory_prefix', help='foo help')
args = parser.parse_args()
print(args)
print(args.file_path)
print(args.directory_prefix)



