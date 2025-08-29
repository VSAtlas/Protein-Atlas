import os
import requests

pdb_ids = [ "6JQR", "4RT7", "5I96", "6ADQ", "6U4J", "6O0K", "5L7I", "1T46", "6WTN", "1IEP", "1OPJ", "3OG7", "5L2I", "4U5J", "2XP2", "4XV2", "1M17", "4AG8", "2GQG", "3CS9", "3OXZ", "5L7D", "3LXK", "4XUF", "2HYY", "4I4E", "3QX3", "2RGC", "4ASD", "3DZY", "5VY4", "3WZE", "3WZD", "4AGC", "4TWP", "5MO4", "4R7H", "6GQO", "5C7X", "2ITN", "3ZOS", "2ZGQ", "2WGJ", "3G0E", "4U2P", "6O0L", "5VA0", "3ZBF", "2E2B", "5TQH"
            #"1IEP", "5MO4", "4R7H", "3WZE", "3WZD", "4AGC", "4TWP", "6GQO", "3ERT", "5C7X"
#
# "4GT3","1OHR", "1HVR", "3CLN", "2P16", "1B9V", "1A30", "1J1Z", "2P54", "3K8Y", "4AKE", "1A6Z", "1GPK"
# ,"6T5U", "1J1Z", "6BD4", "1LAF",
#      "1XCA", "1HL4", "2PRH", "1EYB", "6BD4", "8D0J", "8GBL", "7CR3", "8IJ3", "9J8Z",
#      "7E4T", "8BDC", "7QDS", "9F9Q", "7DF7", "2RXN", "1INS", "4W51", "5H9A", "1GUD",
#      "4W58", "2KMX", "10GS", "4L4R", "9EA7", "8EPH", "1BLR", "2PRL", "2PRM", "7PQT",
#      "2HAU", "3L5G", "5FVT", "6C5R", "4U3Z", "6FYT", "5Z3J", "6NMU", "5LXO", "3C4S",
#      "6L2P", "3O5Y", "5HDQ", "4L2Y", "6Y1S", "4B4E", "2HYY", "5IJO", "3U55", "5Q2W",
#      "6GQ6", "3F8C", "4J7J", "5C0W", "3EMK", "2Z1Y", "6OLQ", "4TNE", "6RXP", "3O2F",
#      "4UYO", "3OZ3", "4CPV", "5T2U", "4W1S", "3KZ9", "5M45", "6A84", "3C87", "6N9T",
#      "5AYF", "6EQU", "4JTI", "6V6G", "3WJ7", "5U57", "6VZI", "3J3K", "5X6F", "4RC6",
#      "6L5R", "5KTO", "6W6X", "3UQV", "5B77", "6J03", "4C2M", "5I6B", "3RTF", "4WHX",
#      "6ZDW", "5ZXZ", "3ZRS", "6Q0B", "4KKF", "5OZ5", "6XKD", "4ZT2", "5MMO", "6UCI",
#      "4L7S", "5OYF", "6Z6Z", "3PXC", "5AGY", "6T6V", "4GZV", "5C7X", "6Y9P"
]




os.makedirs("input_pdbs", exist_ok=True)

for pdb_id in pdb_ids:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url)
    if r.status_code == 200:
        with open(f"input_pdbs/{pdb_id}.pdb", "wb") as f:
            f.write(r.content)
        print(f"Downloaded {pdb_id}")
    else:
        print(f"Failed to download {pdb_id}")
