"""Load private credentials at launch, keeping keys out of MCP config and argv."""
import argparse
import json
import os
from pathlib import Path
import subprocess


def main():
    p = argparse.ArgumentParser()
    p.add_argument('service', choices=['amass', 'nebius', 'japanfold'])
    p.add_argument('--env-file', type=Path)
    p.add_argument('--credentials', type=Path)
    args = p.parse_args()
    env = os.environ.copy()
    if args.service in ['amass', 'japanfold']:
        values = {}
        if args.env_file:
            for line in args.env_file.read_text().splitlines():
                if '=' in line and not line.lstrip().startswith('#'):
                    key, value = line.split('=', 1)
                    values[key.strip()] = value.strip().strip('\"\'')
        name = 'AMASS_API_KEY' if args.service == 'amass' else 'JAPANFOLD_API_KEY'
        key = env.get(name) or values.get(name)
        url = 'https://api.amass.tech/mcp' if args.service == 'amass' else 'https://api.japanfold.aiand.com/mcp'
    else:
        doc = json.loads(args.credentials.read_text())
        api = doc['scientific_api']
        key, url = api['api_key'], api['mcp_url']
    if not key or not url.startswith('https://'):
        raise SystemExit('Missing credential or invalid MCP URL')
    env['HACKATHON_MCP_AUTH'] = 'Bearer ' + key
    # Literal variable reference is expanded by mcp-remote, not the shell.
    cmd = ['npx', '-y', 'mcp-remote@0.16.6', url, '--header',
           'Authorization:${HACKATHON_MCP_AUTH}', '--transport', 'http-only']
    raise SystemExit(subprocess.call(cmd, env=env))

if __name__ == '__main__':
    main()
