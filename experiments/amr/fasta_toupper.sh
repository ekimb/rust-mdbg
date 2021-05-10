 awk 'BEGIN{ cmd = "tr a-z A-Z" } NR % 2 { print } NR % 2 == 0 { print | cmd; close(cmd) }' a

