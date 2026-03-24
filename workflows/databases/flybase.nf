process flybase {
  when: { params.databases.flybase.run }

  output:
  path('*.csv')

  """
  resolve_remote() {
    local remote="${params.databases.flybase.remote}"
    if [[ "\$remote" != *"*"* ]]; then
      printf '%s\n' "\$remote"
      return
    fi

    local index_url="\$(dirname "\$remote")/index.html"
    local match
    match=\$(curl -fsSL "\$index_url" | grep -oE 'ncRNA[^"[:space:]]*\\.json\\.gz' | head -n 1)
    [[ -n "\$match" ]]
    printf '%s/%s\n' "\$(dirname "\$remote")" "\$match"
  }

  wget -O - "\$(resolve_remote)" | gzip -d > flybase.json
  rnac flybase parse flybase.json .
  """
}
