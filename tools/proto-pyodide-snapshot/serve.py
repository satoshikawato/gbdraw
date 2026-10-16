"""PROTOTYPE (proto/pyodide-snapshot, CI measurement only): test server.

Serves <root> at / and <snapshot dir> at /__snapshot__/ (GET; PUT for the
browser-made snapshot), so the snapshot is same-origin with the app but never
inside the repository. With PROTO_SNAPSHOT_URL set, gbdraw/web/js/config.js is
served with PYODIDE_SNAPSHOT_URL set to it (a build that ships the snapshot).
Usage: python3 serve.py <root> <snapshot dir> <port>
"""
import os
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

ROOT, SNAP_DIR, PORT = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), int(sys.argv[3])
PREFIX = '/__snapshot__/'
CONFIG_URL = os.environ.get('PROTO_SNAPSHOT_URL')


class Handler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def _snap_path(self):
        rel = self.path.split('?', 1)[0][len(PREFIX):]
        path = os.path.normpath(os.path.join(SNAP_DIR, rel))
        if not path.startswith(SNAP_DIR + os.sep):
            raise ValueError(rel)
        return path

    def do_GET(self):
        clean = self.path.split('?', 1)[0]
        if clean.startswith(PREFIX):
            try:
                path = self._snap_path()
                with open(path, 'rb') as fh:
                    data = fh.read()
            except (OSError, ValueError):
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header('Content-Type', 'application/octet-stream')
            self.send_header('Content-Length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        if CONFIG_URL and clean == '/gbdraw/web/js/config.js':
            with open(os.path.join(ROOT, 'gbdraw/web/js/config.js'), encoding='utf-8') as fh:
                text = fh.read()
            old = 'export const PYODIDE_SNAPSHOT_URL = null;'
            assert old in text
            data = text.replace(old, f'export const PYODIDE_SNAPSHOT_URL = "{CONFIG_URL}";').encode()
            self.send_response(200)
            self.send_header('Content-Type', 'text/javascript')
            self.send_header('Content-Length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        super().do_GET()

    def do_PUT(self):
        if not self.path.startswith(PREFIX):
            self.send_error(405)
            return
        # PUT <name>?offset=N writes at N (chunked uploads keep each request
        # small enough for the Playwright driver); ?done=1 renames .part.
        path = self._snap_path()
        query = dict(q.split('=', 1) for q in (self.path.split('?', 1)[1].split('&') if '?' in self.path else []))
        length = int(self.headers['Content-Length'])
        os.makedirs(os.path.dirname(path), exist_ok=True)
        offset = int(query.get('offset', '0'))
        mode = 'r+b' if offset and os.path.exists(path + '.part') else 'wb'
        remaining = length
        with open(path + '.part', mode) as fh:
            fh.seek(offset)
            while remaining:
                chunk = self.rfile.read(min(remaining, 1 << 20))
                if not chunk:
                    break
                fh.write(chunk)
                remaining -= len(chunk)
        if query.get('done') == '1':
            os.replace(path + '.part', path)
        self.send_response(201)
        self.send_header('Content-Length', '0')
        self.end_headers()


ThreadingHTTPServer(('127.0.0.1', PORT), partial(Handler, directory=ROOT)).serve_forever()
