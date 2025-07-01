package lib

import groovy.transform.CompileStatic

@CompileStatic
class FastqScan {
    static List<Tuple2<Map<String, String>, String>> findFastqDirs(String rootDir) {
        def fastqDirs = []
        def seenNames = new HashSet<String>()
        new File(rootDir).eachFileRecurse { file ->
            if (file.isDirectory()) {
                def hasFastq = file.listFiles()?.any { f ->
                    f.name ==~ /.*\.(fastq|fq)(\.gz)?$/
                }
                if (hasFastq) {
                    def name = file.name
                    def path = file.absolutePath
                    if (seenNames.contains(name)) {
                        throw new RuntimeException("Name clash: multiple directories named '${name}' found.")
                    }
                    seenNames << name
                    def meta = [barcode: name] as Map<String, String>
                    fastqDirs << new Tuple2(meta, path)
                }
            }
        }

        if (fastqDirs.isEmpty()) {
            throw new RuntimeException("ERROR: No subdirectories with .fastq, .fq, .fastq.gz, or .fq.gz files found in '${rootDir}'")
        }

        return fastqDirs
    }
}
