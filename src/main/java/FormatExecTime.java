import org.jetbrains.annotations.*;
final class FormatExecTime {

    private FormatExecTime() {
        throw new AssertionError("Utility class should not be instantiated!");
    }

    /**
     * Converts an elapsed time in nanoseconds to a human-readable string
     * using the largest appropriate time unit.
     *
     * @param ns elapsed time in nanoseconds
     * @return a string like "42 milliseconds", "3.5 seconds", or "0.2 minutes"
     */
    public static @NotNull String formatDuration(long ns) {
        if (ns >= 60_000_000_000L) {
            double minutes = ns / 60_000_000_000.0;
            return String.format("%.3f minutes", minutes);
        } else if (ns >= 1_000_000_000L) {
            double seconds = ns / 1_000_000_000.0;
            return String.format("%.3f seconds", seconds);
        } else if (ns >= 1_000_000L) {
            double ms = ns / 1_000_000.0;
            return String.format("%.3f milliseconds", ms);
        } else if (ns >= 1_000L) {
            double us = ns / 1_000.0;
            return String.format("%.3f microseconds", us);
        } else {
            return ns + " nanoseconds";
        }
    }
}